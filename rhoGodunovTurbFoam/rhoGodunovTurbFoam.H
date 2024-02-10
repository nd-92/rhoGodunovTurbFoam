#include "readTimeControls.H"
Co.calculate();
Co.setDeltaT(maxDeltaT, maxCo, adjustTimeStep);

runTime++;
Info << "\nTime = " << runTime.value() << "\n";

if (useOldApproach)
{
    // Solve the approximate Riemann problem for this time step
    dbnsFlux.computeFlux();

    // Time integration
    solve(                                                                                                                                                                                            //
        fvm::ddt(rho)                                                                                                                                                                                 //
        + fvc::div(dbnsFlux.rhoFlux()));                                                                                                                                                              //
    solve(                                                                                                                                                                                            //
        fvm::ddt(rhoU)                                                                                                                                                                                //
        + fvc::div(dbnsFlux.rhoUFlux())                                                                                                                                                               //
        - fvc::laplacian(turbulence->muEff(), U)                                                                                                                                                      //
        - fvc::div(turbulence->muEff() * dev2(Foam::T(fvc::grad(U)))));                                                                                                                               //
    solve(                                                                                                                                                                                            //
        fvm::ddt(rhoE)                                                                                                                                                                                //
        + fvc::div(dbnsFlux.rhoEFlux())                                                                                                                                                               //
        - fvc::div((fvc::interpolate(turbulence->muEff()) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), turbulence->muEff() * dev2(Foam::T(fvc::grad(U))))) & fvc::interpolate(U)) //
        - fvc::laplacian(turbulence->alphaEff(), e));                                                                                                                                                 //

    // Compute U
    U = rhoU / rho;
    U.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

    // Calculate internal energy
    e = ((rhoE / rho) - (0.5 * magSqr(U)));
    e.correctBoundaryConditions();

    // Correct thermodynamics
    thermo.correct();
    rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

    // Compute p from rho
    thermo.p().ref() = rho() / thermo.psi();
    thermo.p().correctBoundaryConditions();
    rho.boundaryFieldRef() == thermo.psi().boundaryField() * thermo.p().boundaryField();

    // Correct turbulence fields
    turbulence->correct();

    // Update the shock sensor
    shockDetector.update();
}
else
{
    // Solve the approximate Riemann problem for this time step
    dbnsFlux.computeFlux();

    // Explicit integration of rho
    solve(fvm::ddt(rho) + fvc::div(dbnsFlux.rhoFlux()));

    // Explicit integration of rhoU
    solve(fvm::ddt(rhoU) + fvc::div(dbnsFlux.rhoUFlux()));

    // Assign U
    U.ref() = rhoU() / rho();
    U.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

    // Account for viscous terms in momentum
    solve(fvm::ddt(rho, U) - fvc::ddt(rho, U) - fvm::laplacian(turbulence->muEff(), U) - fvc::div(turbulence->muEff() * dev2(Foam::T(fvc::grad(U)))));
    rhoU = rho * U;

    // Solve energy equation
    solve(fvm::ddt(rhoE) + fvc::div(dbnsFlux.rhoEFlux()) - fvc::div((fvc::interpolate(turbulence->muEff()) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), turbulence->muEff() * dev2(Foam::T(fvc::grad(U))))) & fvc::interpolate(U)));

    // Assign e
    e = (rhoE / rho) - (0.5 * magSqr(U));
    e.correctBoundaryConditions();
    thermo.correct();
    rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + (0.5 * magSqr(U.boundaryField())));

    // Account for viscous terms in energy
    solve(fvm::ddt(rho, e) - fvc::ddt(rho, e) - fvm::laplacian(turbulence->alphaEff(), e));

    // Update thermo properties
    thermo.correct();

    // Finish time step
    rhoE = rho * (e + (0.5 * magSqr(U)));
    thermo.p().ref() = rho() / thermo.psi();
    thermo.p().correctBoundaryConditions();
    rho.boundaryFieldRef() == thermo.psi().boundaryField() * thermo.p().boundaryField();

    // Update turbulence fields
    turbulence->correct();

    // Update the shock sensor
    shockDetector.update();
}
