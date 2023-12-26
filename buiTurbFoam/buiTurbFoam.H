#include "readTimeControls.H"
Co.calculate();
Co.setDeltaT(maxDeltaT, maxCo, adjustTimeStep);

runTime++;
Info << "\nTime = " << runTime.value() << "\n";

// Low storage Runge-Kutta time integration
rkLoop(i)
{
    // Solve the approximate Riemann problem for this time step
    dbnsFlux.computeFlux();

    // Assemble and solve the density equation
    solve(                                   //
        1 / rkCoeffs.beta(i) * fvm::ddt(rho) //
        + fvc::div(dbnsFlux.rhoFlux()));     //

    // Assemble and solve the momentum equation
    solve(                                                              //
        1 / rkCoeffs.beta(i) * fvm::ddt(rhoU)                           //
        + fvc::div(dbnsFlux.rhoUFlux())                                 //
        - fvc::laplacian(turbulence->muEff(), U)                        //
        - fvc::div(turbulence->muEff() * dev2(Foam::T(fvc::grad(U))))); //

    // Assemble and solve the energy equation
    solve(                                                                                                                                                                                            //
        1 / rkCoeffs.beta(i) * fvm::ddt(rhoE)                                                                                                                                                         //
        + fvc::div(dbnsFlux.rhoEFlux())                                                                                                                                                               //
        - fvc::div((fvc::interpolate(turbulence->muEff()) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), turbulence->muEff() * dev2(Foam::T(fvc::grad(U))))) & fvc::interpolate(U)) //
        - fvc::laplacian(turbulence->alphaEff(), e));                                                                                                                                                 //

    // Compute U
    U = rhoU / rho;
    U.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

    // Calculate internal energy
    // e = (rhoE / rho) - (0.5 * magSqr(U));
    e = (rhoE / rho) - (0.5 * magSqr(U)) - turbulence->k();
    e.correctBoundaryConditions();

    // Correct thermodynamics
    thermo.correct();
    rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5 * magSqr(U.boundaryField()) + turbulence->k().boundaryField());

    // Compute p from rho
    p.ref() = rho() / thermo.psi();
    p.correctBoundaryConditions();
    rho.boundaryFieldRef() == thermo.psi().boundaryField() * p.boundaryField();
}

// Correct turbulence fields
turbulence->correct();
