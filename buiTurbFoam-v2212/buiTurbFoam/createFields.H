Info << "Reading thermophysical properties\n"
     << endl;
autoPtr<psiThermo> pThermo(psiThermo::New(mesh));
psiThermo &thermo = pThermo();
thermo.validate(args.executable(), "h");
volScalarField &he = thermo.he();
volScalarField &p = thermo.p();
volScalarField &T = thermo.T();

volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.rho());

Info << "Reading fields U, rho * U, rho * E\n"
     << endl;

volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
volVectorField rhoU(
    IOobject(
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    rho *U);
volScalarField rhoE(
    IOobject(
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    rho *(he + 0.5 * magSqr(U)));

Info << "Reading field upwindingFactor" << endl;
const volScalarField upwindingFactor(
    IOobject(
        "upwindingFactor",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(
        "one",
        dimless,
        1));
const surfaceScalarField buiEps(
    IOobject(
        "buiEps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    fvc::interpolate(upwindingFactor));

autoPtr<basicNumericFlux> dbnsFluxPtr = basicNumericFlux::New(
    p,
    U,
    T,
    buiEps,
    thermo);
basicNumericFlux &dbnsFlux = dbnsFluxPtr();

const surfaceScalarField &phi = dbnsFlux.rhoFlux();

Info << "Creating turbulence model" << endl;
autoPtr<compressible::turbulenceModel> turbulence(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        thermo));

volScalarField muEff(
    "muEff",
    turbulence->muEff());

volTensorField tauMC(
    "tauMC",
    muEff *dev2(Foam::T(fvc::grad(U))));

surfaceScalarField sigmaDotU(
    "sigmaDotU",
    (fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), tauMC)) & fvc::interpolate(U));

// Runge-Kutta coefficient
scalarList beta(4);
beta[0] = 0.1100;
beta[1] = 0.2766;
beta[2] = 0.5000;
beta[3] = 1.0000;

// Constants for blending
const dimensionedScalar T_inf(
    "T_inf",
    T.dimensions(),
    readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("T")));
const dimensionedScalar p_inf(
    "p_inf",
    p.dimensions(),
    readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("p")));
const dimensionedScalar rho_inf(
    "rho_inf",
    rho.dimensions(),
    readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("rho")));
const dimensionedVector U_inf(
    "U_inf",
    U.dimensions(),
    vector(mesh.solutionDict().subDict("freestreamProperties").lookup("U")));

// Blending coefficient
const volScalarField acousticBlending(
    IOobject(
        "acousticBlending",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(
        "one",
        dimless,
        1));

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
