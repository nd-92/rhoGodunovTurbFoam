{
    // Compute U
    U = rhoU / rho;
    U.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

    // Calculate enthalpy from rhoE
    e.ref() = (rhoE / rho) - (0.5 * magSqr(U));
    e.correctBoundaryConditions();

    // Correct thermodynamics
    thermo.correct();
    rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

    // Compute p from rho
    p.ref() = rho() / thermo.psi();
    p.correctBoundaryConditions();
    rho.boundaryFieldRef() == thermo.psi().boundaryField() * p.boundaryField();
}
