// Thermo fields
Info << "Reading thermophysical properties" << endl;
autoPtr<psiThermo> pThermo(psiThermo::New(mesh));
psiThermo &thermo = pThermo();
// thermo.validate(args.executable(), "h");
volScalarField &e = thermo.he();
volScalarField &p = thermo.p();
volScalarField &T = thermo.T();

// Momentum fields
Info << "Reading fields rho, U, rho * U, rho * E" << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    thermo.rho());
volVectorField rhoU(
    IOobject(
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.rho() * U);
volScalarField rhoE(
    IOobject(
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.rho() * (e + 0.5 * magSqr(U)));

// Upwinding factor
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

// Numeric flux
autoPtr<basicNumericFlux> dbnsFluxPtr = basicNumericFlux::New(
    p,
    U,
    T,
    upwindingFactor,
    thermo);
basicNumericFlux &dbnsFlux = dbnsFluxPtr();
const surfaceScalarField &phi = dbnsFlux.rhoFlux();

// Turbulence fields
Info << "Creating turbulence model" << endl;
autoPtr<compressible::turbulenceModel> turbulence(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        thermo));
