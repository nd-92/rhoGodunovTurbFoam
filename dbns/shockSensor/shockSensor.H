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

Class
    shockSensor

Description
    Implementation of a density gradient based shock sensor

SourceFiles
    shockSensor.H

\*---------------------------------------------------------------------------*/

#ifndef shockSensor_H
#define shockSensor_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Helper function to construct local shock sensor DT
    inline const volScalarField localFilterDT(
        // inline const dimensionedScalar localFilterDT(
        const fvMesh &mesh,
        const Time &runTime)
    {
        const scalar filterScale = mesh.schemesDict().subDict("divSchemes").subDict("dbns").lookupOrDefault("shockFilter", 1.0);

        volScalarField filterWidth(
            IOobject(
                "shockSensorFilterWidth",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("zero", dimLength, 0));

        // Loop over all cells, setting volume
        forAll(mesh.cells(), i)
        {
            filterWidth[i] = pow(mesh.V()[i], 1.0 / 3.0) * filterScale;
        }

        // Loop over all boundary faces, setting area
        volScalarField::Boundary &filterWidthBf = filterWidth.boundaryFieldRef();
        forAll(filterWidthBf, patchi)
        {
            // Test setting area
            filterWidthBf[patchi] = sqrt(mesh.magSf().boundaryField()[patchi]) * filterScale;
        }

        return sqr(filterWidth) / dimensionedScalar("dt_", dimTime, mesh.time().deltaTValue());
    };

    // Helper function to construct smoothing switch
    inline bool isSmoothing(
        const fvMesh &mesh)
    {
        const scalar filterScale = mesh.schemesDict().subDict("divSchemes").subDict("dbns").lookupOrDefault("shockFilter", 1.0);
        if (filterScale > 0)
        {
            Info << "Initialised shock sensor with local filter with scale = " << filterScale << " : 1" << endl;
            return true;
        }
        else
        {
            return false;
        }
    };

    /*---------------------------------------------------------------------------*\
                               Class shockSensor Declaration
    \*---------------------------------------------------------------------------*/

    template <class Thermo>
    class shockSensor
    {
    public:
        // Constructor
        inline shockSensor(
            const volVectorField &U,
            const fvMesh &mesh,
            const Thermo &thermo,
            const Time &runTime)
            : U_(U),
              mesh_(mesh),
              thermo_(thermo),
              runTime_(runTime),
              DT_(localFilterDT(mesh_, runTime_)),
              shockSensor_(
                  IOobject(
                      "buiEps",
                      runTime.timeName(),
                      mesh,
                      IOobject::MUST_READ,
                      IOobject::AUTO_WRITE),
                  mesh){};

        // Destructor
        inline ~shockSensor(){};

        // Update the shock sensor
        inline void update()
        {
            // Shock sensor is proportional to mag(grad(rho))
            shockSensor_ = mag(fvc::grad(thermo_.rho())) / max(mag(fvc::grad(thermo_.rho())));

            // Smoothen the shock sensor by solving the Laplace equation in pseudo-time
            solve(fvm::ddt(shockSensor_) - fvm::laplacian(DT_, shockSensor_));
            shockSensor_.correctBoundaryConditions();

            // Clamp to a range of 0 to 1
            shockSensor_ = (shockSensor_ - min(shockSensor_)) / (max(shockSensor_) - min(shockSensor_));
        }

        // Public const access to the shock sensor
        inline const volScalarField &upwindingFactor() const
        {
            return shockSensor_;
        }

    private:
        // Const access to velocity field
        const volVectorField &U_;

        // Const access to mesh
        const fvMesh &mesh_;

        // Const access to thermo
        const Thermo &thermo_;

        // Const access to runTime
        const Time &runTime_;

        // ddt(Laplacian) operator scale
        const volScalarField DT_;

        // Virgin smoother field
        volScalarField shockSensor_;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //