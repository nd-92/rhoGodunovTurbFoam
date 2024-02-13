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
    Density-based compressible transient solver using the flux difference
    splitting scheme of Bui. Primarily designed for LES.

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
#include "shockSensor.H"
#include "numericFluxes.H"

bool readSolutionApproach(const fvMesh &mesh)
{
    word solutionApproach = "new";

    if (mesh.schemesDict().readIfPresent("solutionApproach", solutionApproach))
    {
        if (solutionApproach == "old")
        {
            Info << "solutionApproach: " << solutionApproach << endl;
            return true;
        }
        else if (solutionApproach == "new")
        {
            Info << "solutionApproach: " << solutionApproach << endl;
            return false;
        }
        else
        {
            FatalErrorInFunction
                << "solutionApproach: " << solutionApproach
                << " is not a valid choice. "
                << "Options are: old, new"
                << abort(FatalError);
            return false;
        }
    }
    return false;
};

int main(int argc, char *argv[])
{

    argList::addNote(
        "Density-based compressible transient solver using the flux difference"
        "splitting scheme of Bui. Primarily designed for LES.");

#include "postProcess.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

#include "createTimeControls.H"

    // const bool useOldApproach = readSolutionApproach(mesh);

    // Get reference values for damping
    const dimensionedScalar T_inf(
        "T_inf",
        thermo.T().dimensions(),
        readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("T")));
    const dimensionedScalar p_inf(
        "p_inf",
        thermo.p().dimensions(),
        readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("p")));
    const dimensionedScalar rho_inf(
        "rho_inf",
        rho.dimensions(),
        readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("rho")));
    const dimensionedVector U_inf(
        "U_inf",
        U.dimensions(),
        vector(mesh.solutionDict().subDict("freestreamProperties").lookup("U")));

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
