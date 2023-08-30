// Compute U
U.ref() = rhoU() / rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

// Bound density
rho = bound(rho, rhoMin);

// Bound temperature
T = bound(T, TMin);

// Calculate enthalpy
he.ref() = rhoE / rho - 0.5 * magSqr(U);
he.correctBoundaryConditions();

// Bound enthalpy
he = bound(he, heMin);

// Correct thermodynamics
thermo.correct();
rhoE.boundaryFieldRef() == rho.boundaryField() * (he.boundaryField() + 0.5 * magSqr(U.boundaryField()));

// Compute p from rho
p.ref() = rho() / thermo.psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == thermo.psi().boundaryField() * p.boundaryField();

// Bound pressure
p = bound(p, pMin);
