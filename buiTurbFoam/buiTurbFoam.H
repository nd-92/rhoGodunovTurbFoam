#include "readTimeControls.H"
#include "compressibleCourantNo.H"
#include "setDeltaT.H"

runTime++;
Info << "\n Time = " << runTime.value() << endl;

// Low storage Runge-Kutta time integration
forAll(beta, i)
{
    // Solve the approximate Riemann problem for this time step
    dbnsFlux.computeFlux();

    // Time integration
    solve(                                                                                                                                                                                            //
        1.0 / beta[i] * fvm::ddt(rho)                                                                                                                                                                 //
        + fvc::div(dbnsFlux.rhoFlux()));                                                                                                                                                              //
    solve(                                                                                                                                                                                            //
        1.0 / beta[i] * fvm::ddt(rhoU)                                                                                                                                                                //
        + fvc::div(dbnsFlux.rhoUFlux())                                                                                                                                                               //
        - fvc::laplacian(turbulence->muEff(), U)                                                                                                                                                      //
        - fvc::div(turbulence->muEff() * dev2(Foam::T(fvc::grad(U)))));                                                                                                                               //
    solve(                                                                                                                                                                                            //
        1.0 / beta[i] * fvm::ddt(rhoE)                                                                                                                                                                //
        + fvc::div(dbnsFlux.rhoEFlux())                                                                                                                                                               //
        - fvc::div((fvc::interpolate(turbulence->muEff()) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), turbulence->muEff() * dev2(Foam::T(fvc::grad(U))))) & fvc::interpolate(U)) //
        - fvc::laplacian(turbulence->alphaEff(), he));                                                                                                                                                //

    // Update fields to new time step
#include "updateFields.H"
}
