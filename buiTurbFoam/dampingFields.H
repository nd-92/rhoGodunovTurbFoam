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
U = (U * acousticBlending) - (U_inf * (acousticBlending - 1));
rho = (rho * acousticBlending) - (rho_inf * (acousticBlending - 1));
rhoU = (rhoU * acousticBlending) - (rho_inf * U_inf * (acousticBlending - 1));
T = (T * acousticBlending) - (T_inf * (acousticBlending - 1));
p = (p * acousticBlending) - (p_inf * (acousticBlending - 1));
