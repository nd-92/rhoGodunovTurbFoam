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
    Density-based compressible transient solver using the flux difference
    splitting scheme of Bui and Runge Kutta 4-stage time integration.
    Primarily designed for LES, and with optional absorption-based acoustic
    damping

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "MDLimiter.H"
#include "firstOrderLimiter.H"
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H"
#include "WangLimiter.H"
#include "numericFlux.H"
#include "buiFlux.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "bound.H"

int main(int argc, char *argv[])
{

    argList::addNote(
        "Density-based compressible transient solver using the flux"
        " difference splitting scheme of Bui and Runge Kutta 4-stage"
        " time integration.  Primarily designed for LES, and with"
        " optional absorption-based acoustic damping.");

#include "postProcess.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

    // Runge-Kutta coefficient
    constexpr double betaReciprocal[4] = {1 / 0.1100, 1 / 0.2766, 2.000, 1.000};

    // Damping switch
    const bool applyDamping = readBool(runTime.controlDict().lookup("applyDamping"));

    // Courant number
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar velMag = 0.0;

#include "createTimeControls.H"

    Info << "Starting time loop" << endl;

    if (applyDamping == true)
    {
#include "dampingFields.H"
        while (runTime.run())
        {
            // Execute main solver loop
#include "buiTurbFoamDamping.H"

            // Write runtime output
            runTime.write();
            runTime.printExecutionTime(Info);
        }

        Info << "End" << endl;

        return 0;
    }
    else
    {
        while (runTime.run())
        {
            // Execute main solver loop
#include "buiTurbFoam.H"

            // Write runtime output
            runTime.write();
            runTime.printExecutionTime(Info);
        }

        Info << "End" << endl;

        return 0;
    }
}

// ************************************************************************* //
