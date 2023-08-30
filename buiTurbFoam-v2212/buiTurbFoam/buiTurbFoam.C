/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buiTurbFoam

Description
    Density-based compressible flow solver based on the Roe flux difference
    splitting scheme of Bui with absorbing sponge zones

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "MDLimiter.H"
#include "firstOrderLimiter.H"
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H"
#include "numericFlux.H"
#include "buiFlux.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "bound.H"

int main(int argc, char *argv[])
{

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
#include "createTimeControls.H"
#include "readFieldBounds.H"

    Info << "\nStarting time loop" << endl;

    while (runTime.run())
    {

#include "readTimeControls.H"
#include "acousticCourantNo.H"
#include "setDeltaT.H"

        runTime++;
        Info << "\n Time = " << runTime.value() << endl;

        // Low storage Runge-Kutta time integration
        forAll(beta, i)
        {
            // Solve the approximate Riemann problem for this time step
            dbnsFlux.computeFlux();

            // Time integration
            solve(1.0 / beta[i] * fvm::ddt(rho) + fvc::div(dbnsFlux.rhoFlux()));
            muEff = turbulence->muEff();
            tauMC = muEff * dev2(Foam::T(fvc::grad(U)));
            solve(1.0 / beta[i] * fvm::ddt(rhoU) + fvc::div(dbnsFlux.rhoUFlux()) - fvc::laplacian(muEff, U) - fvc::div(tauMC));
            sigmaDotU = (fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), tauMC)) & fvc::interpolate(U);
            solve(1.0 / beta[i] * fvm::ddt(rhoE) + fvc::div(dbnsFlux.rhoEFlux()) - fvc::div(sigmaDotU) - fvc::laplacian(turbulence->alphaEff(), he));

            // Update fields to new time step
#include "updateFields.H"
        }

        // Apply acoustic blending
        U = (U * acousticBlending) - (U_inf * (acousticBlending - 1));
        thermo.rho() = (thermo.rho() * acousticBlending) - (rho_inf * (acousticBlending - 1));
        thermo.p() = (thermo.p() * acousticBlending) - (p_inf * (acousticBlending - 1));
        thermo.T() = (thermo.T() * acousticBlending) - (T_inf * (acousticBlending - 1));
        rho = thermo.rho();
        p = thermo.p();
        T = thermo.T();

        // Correct turbulence fields
        turbulence->correct();

        // Write runtime output
        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info << "End \n";

    return 0;
}

// ************************************************************************* //
