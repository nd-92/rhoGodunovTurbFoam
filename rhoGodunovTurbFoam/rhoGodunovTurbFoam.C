/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    rhoGodunovTurbFoam

Description
    Implicit density-based compressible transient solver using the flux
    difference splitting scheme of Bui. Primarily designed for LES.

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
#include "rungeKutta.H"
#include "acousticCourantNo.H"

int main(int argc, char *argv[])
{

    argList::addNote(
        "Implicit density-based compressible transient solver using the flux"
        "difference splitting scheme of Bui. Primarily designed for LES.");

#include "postProcess.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

    // Runge-Kutta coefficient
    // rungeKutta rkCoeffs;

    // Acoustic Courant number CFL control
    acousticCourantNo<psiThermo> Co(U, mesh, thermo, runTime);

#include "createTimeControls.H"

    Info << "Starting time loop" << endl;

    while (runTime.run())
    {
        // Execute main solver loop
#include "rhoGodunovTurbFoam.H"

        // Write runtime output
        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info << "End" << endl;

    return 0;
}

// ************************************************************************* //