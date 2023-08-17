// Compute U
//    U = rhoU/rho;
U.ref() = rhoU() / rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

// Calculate enthalpy from rhoE
//    const volScalarField Cp = thermo.Cp();
//    const volScalarField Cv = thermo.Cv();

//     h = rhoE/rho - 0.5*magSqr(U) + p/rho;
// Alternative formulation, Felipe Alves Portela TU Delft
// h.ref() = thermo.Cp()/thermo.Cv()*(rhoE/rho - 0.5*magSqr(U));
// h.correctBoundaryConditions();
e.ref() = rhoE / rho - 0.5 * magSqr(U);
e.correctBoundaryConditions();

// Correct thermodynamics
thermo.correct();
rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

// Compute p from rho
//    p = rho*(thermo.Cp() - thermo.Cv())*T;
p.ref() = rho() / thermo.psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == thermo.psi().boundaryField() * p.boundaryField();
