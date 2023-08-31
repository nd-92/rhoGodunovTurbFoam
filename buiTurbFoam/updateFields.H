// // Compute U
// U.ref() = rhoU() / rho();
// U.correctBoundaryConditions();
// rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

// // Bound density
// rho = bound(rho, rhoMin);

// // Bound temperature
// T = bound(T, TMin);

// // Calculate enthalpy
// he.ref() = rhoE / rho - 0.5 * magSqr(U);
// he.correctBoundaryConditions();

// // Bound enthalpy
// // he = bound(he, heMin);

// // Correct thermodynamics
// thermo.correct();
// rhoE.boundaryFieldRef() == rho.boundaryField() * (he.boundaryField() + 0.5 * magSqr(U.boundaryField()));

// // Compute p from rho
// p.ref() = rho() / thermo.psi();
// p.correctBoundaryConditions();
// rho.boundaryFieldRef() == thermo.psi().boundaryField() * p.boundaryField();

// // Bound pressure
// // p = bound(p, pMin);

{
    // Compute U
    U = rhoU / rho;
    U.correctBoundaryConditions();

    // Calculate enthalpy from rhoE
    // const volScalarField Cp = thermo->Cp();
    // const volScalarField Cv = thermo->Cv();

    // h = rhoE/rho - 0.5*magSqr(U) + p/rho;
    // Alternative formulation, Felipe Alves Portela TU Delft
    he = thermo.Cp() / thermo.Cv() * (rhoE / rho - 0.5 * magSqr(U));
    he.correctBoundaryConditions();

    // Bound enthalpy
    // dimensionedScalar CpMin = min(Cp);
    // dimensionedScalar CpMax = max(Cp);

    // dimensionedScalar hMin = CpMin * TMin;
    // dimensionedScalar hMax = CpMax * TMax;

    he = bound(he, heMin);

    // Correct thermodynamics
    thermo.correct();

    // Bound density
    rho = bound(rho, rhoMin);

    // Compute p from rho
    p = rho * (thermo.Cp() - thermo.Cv()) * T;
    p.correctBoundaryConditions();

    // Bound pressure
    p = bound(p, pMin);
}
