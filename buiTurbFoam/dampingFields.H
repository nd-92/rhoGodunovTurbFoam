// Get reference values for damping
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
const dimensionedScalar e_inf(
    "e_inf",
    average(thermo.he(p_inf + (0 * thermo.p()), T_inf + (0 * thermo.T()))));

// Get freestream effective viscosity and thermal diffusivity
const bool turbulenceDamping = readBool(runTime.controlDict().lookup("turbulenceDamping"));
// if (turbulenceDamping == true)
// {
const dimensionedScalar muEff_inf(
    "muEff_inf",
    turbulence->muEff()->dimensions(),
    readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("muEff")));
const dimensionedScalar alphaEff_inf(
    "alphaEff_inf",
    turbulence->alphaEff()->dimensions(),
    readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("alphaEff")));
//}

// Read blending coefficient, optionally assign = 1
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

// Initially damp fields
rho = (rho * acousticBlending) - (rho_inf * (acousticBlending - 1));
rhoU = (rhoU * acousticBlending) - (rho_inf * U_inf * (acousticBlending - 1));
rhoE = (rhoE * acousticBlending) - (rho_inf * e_inf * (acousticBlending - 1));
U = (U * acousticBlending) - (U_inf * (acousticBlending - 1));
e = (e * acousticBlending) - (e_inf * (acousticBlending - 1));
T = (T * acousticBlending) - (T_inf * (acousticBlending - 1));
p = (p * acousticBlending) - (p_inf * (acousticBlending - 1));
if (turbulenceDamping == true)
{
    turbulence->muEff() = (turbulence->muEff() * acousticBlending) - (muEff_inf * (acousticBlending - 1));
    turbulence->alphaEff() = (turbulence->alphaEff() * acousticBlending) - (alphaEff_inf * (acousticBlending - 1));
}
