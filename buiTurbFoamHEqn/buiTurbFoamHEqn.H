#include "readTimeControls.H"
#include "acousticCourantNo.H"
#include "setDeltaT.H"

runTime++;
Info << "\nTime = " << runTime.value() << "\n";

// Low storage Runge-Kutta time integration
// forAll(beta, i)
for (uint8_t i = 0; i < 4; i++)
{
    // Solve the approximate Riemann problem for this time step
    dbnsFlux.computeFlux();

    // Time integration
    solve(                                   //
        1 / rkCoeffs.beta(i) * fvm::ddt(rho) //
        + fvc::div(dbnsFlux.rhoFlux()));     //

    solve(                                     //
        1 / rkCoeffs.beta(i) * fvm::ddt(rhoU)  //
        + fvc::div(dbnsFlux.rhoUFlux())        //
        + fvc::div(turbulence->devRhoReff())); //

    solve(
        1.0 / rkCoeffs.beta(i) * fvm::ddt(rhoE)       //
        + fvc::div(dbnsFlux.rhoEFlux())               //
        + fvc::div(turbulence->devRhoReff() & U)      //
        - fvc::laplacian(turbulence->alphaEff(), h)); //

    // Compute U
    U = rhoU / rho;
    U.correctBoundaryConditions();

    // Calculate enthalpy
    h = thermo.gamma() * ((rhoE / rho) - (0.5 * magSqr(U)));
    h.correctBoundaryConditions();

    // Correct thermodynamics
    thermo.correct();

    // Compute p from rho
    p = rho * (thermo.Cp() - thermo.Cv()) * thermo.T();
    p.correctBoundaryConditions();

    // Correct turbulence fields
    turbulence->correct();
}
