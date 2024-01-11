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
    acousticCourantNo

Description
    Implementation of the acoustic Courant number

SourceFiles
    acousticCourantNo.H

\*---------------------------------------------------------------------------*/

#ifndef acousticCourantNo_H
#define acousticCourantNo_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Helper function to return minimum cell dimension
    // inline const scalarField minDeltaxyz(
    //     const fvMesh &mesh,
    //     const Time &runTime)
    // {

    //     // Initialise the volScalarField
    //     // volScalarField cellDim(
    //     //     IOobject(
    //     //         "cellDim",
    //     //         runTime.timeName(),
    //     //         mesh,
    //     //         IOobject::NO_READ,
    //     //         IOobject::NO_WRITE),
    //     //     mesh,
    //     //     dimensionedScalar(
    //     //         "one",
    //     //         dimLength,
    //     //         1));

    //     const cellList &cells = mesh.cells();
    //     const vectorField &cellC = mesh.cellCentres();
    //     const vectorField &faceC = mesh.faceCentres();
    //     const vectorField faceN(mesh.faceAreas() / mag(mesh.faceAreas()));

    //     scalarField cellDim(mesh.nCells());

    //     // Loop over all cells
    //     forAll(cells, celli)
    //     {
    //         // Initial guess of distance from cell centroid to face 0
    //         scalar deltaMaxTmp = mag(faceN[0] & (faceC[0] - cellC[celli]));
    //         const labelList &cFaces = cells[celli];
    //         const point &cc = cellC[celli];

    //         // Loop over faces of cell
    //         forAll(cFaces, cFacei)
    //         {
    //             const label facei = cFaces[cFacei];
    //             const point &fc = faceC[facei];
    //             const vector &n = faceN[facei];

    //             // Distance from face to cell centroids
    //             const scalar tmp = mag(n & (fc - cc));

    //             // Assign min value
    //             if (tmp < deltaMaxTmp)
    //             {
    //                 deltaMaxTmp = tmp;
    //             }
    //         }

    //         // Cell dimension is 2 * distance from cell centroid to face centroid
    //         cellDim[celli] = 2 * deltaMaxTmp;
    //     }

    //     return cellDim;
    // };

    /*---------------------------------------------------------------------------*\
                            Class acousticCourantNo Declaration
    \*---------------------------------------------------------------------------*/

    // template <class Thermo>
    // class acousticCourantNo
    // {
    // public:
    //     // Constructor
    //     inline acousticCourantNo(
    //         const volVectorField &U,
    //         const fvMesh &mesh,
    //         const Thermo &thermo,
    //         Time &runTime)
    //         : U_(U),
    //           mesh_(mesh),
    //           thermo_(thermo),
    //           runTime_(runTime),
    //           cellDim_(minDeltaxyz(mesh, runTime)),
    //           fc_(0 * minDeltaxyz(mesh, runTime)),
    //           nCells_(static_cast<scalar>(mesh_.nCells())),
    //           CoNum_(0),
    //           meanCoNum_(0),
    //           gamma_(thermo.gamma()),
    //           Cp_(thermo.Cp()),
    //           Cv_(thermo.Cv())
    //     {
    //         // Only need to evaluate this once at construction time
    //         // Avoid evaluation every time update is called
    //         if (mesh_.nInternalFaces() == false)
    //         {
    //             throw std::invalid_argument("No internal mesh faces.");
    //         }
    //     };

    //     // Destructor
    //     inline ~acousticCourantNo(){};

    //     // Calculate and print acoustic Courant number
    //     inline void calculate()
    //     {
    //         // Calculate dimensionless acoustic frequency
    //         // fc_ = (mag(U_.primitiveField()) + sqrt(gamma_.primitiveField() * (Cp_.primitiveField() - Cv_.primitiveField()) * thermo_.T().primitiveField())) / cellDim_;

    //         // CoNum_ = 0;
    //         // meanCoNum_ = 0;
    //         // forAll(mesh_.cells(), i)
    //         // {
    //         //     meanCoNum_ = meanCoNum_ + fc_[i];
    //         //     if (fc_[i] > CoNum_)
    //         //     {
    //         //         CoNum_ = fc_[i];
    //         //     }
    //         // }
    //         // CoNum_ = CoNum_ * runTime_.deltaTValue();
    //         // meanCoNum_ = meanCoNum_ / nCells_;

    //         // Info << max(fc_) << endl;

    //         // Get max Courant number
    //         // const dimensionedScalar fcMax = gMax(fc_);
    //         // CoNum_ = gMax(fc_.primitiveField()) * runTime_.deltaTValue();

    //         // Get mean Courant number
    //         // const dimensionedScalar fcMean = gSum(fc_);
    //         // meanCoNum_ = gSum(fc_.primitiveField()) * runTime_.deltaTValue() / nCells_;

    //         // Output
    //         // Info << "Mean and max acoustic Courant Numbers = " << meanCoNum_ << " " << CoNum_ << endl;
    //     };

    //     // Update the time step
    //     inline void setDeltaT(
    //         const scalar maxDeltaT,
    //         const scalar maxCo,
    //         const bool adjustTimeStep)
    //     {
    //         if (adjustTimeStep)
    //         {
    //             const scalar maxDeltaTFact = maxCo / (CoNum_ + SMALL);
    //             const scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2);

    //             runTime_.setDeltaT(
    //                 min(
    //                     deltaTFact * runTime_.deltaTValue(),
    //                     maxDeltaT));

    //             Info << "deltaT = " << runTime_.deltaTValue() << endl;
    //         }
    //     }

    // private:
    //     // Const access to velocity field
    //     const volVectorField &U_;

    //     // Const access to mesh
    //     const fvMesh &mesh_;

    //     // Const access to thermo
    //     const Thermo &thermo_;

    //     // Non-const access to runTime
    //     Time &runTime_;

    //     // Mesh cell dimension
    //     const scalarField cellDim_;

    //     // Cell acoustic frequency
    //     scalarField fc_;

    //     // Scalar number of mesh cells
    //     const scalar nCells_;

    //     // Max Courant number
    //     scalar CoNum_;

    //     // Mean Courant number
    //     scalar meanCoNum_;

    //     // Ratio of specific heats
    //     const volScalarField &gamma_;

    //     // Specific heat capacities
    //     const volScalarField &Cp_;
    //     const volScalarField &Cv_;
    // };

    /*---------------------------------------------------------------------------*\
                            Class acousticCourantNo Declaration
    \*---------------------------------------------------------------------------*/

    template <class Thermo>
    class acousticCourantNo
    {
    public:
        // Constructor
        inline acousticCourantNo(
            const volVectorField &U,
            const fvMesh &mesh,
            const Thermo &thermo,
            Time &runTime)
            : U_(U),
              mesh_(mesh),
              thermo_(thermo),
              runTime_(runTime),
              amaxSf_(0.5 * (mag((fvc::interpolate(U_) & mesh_.Sf())) + (fvc::interpolate(sqrt(thermo_.gamma() * (thermo_.Cp() - thermo_.Cv()) * thermo_.T())) * mesh_.magSf()))),
              sumAmaxSf_(fvc::surfaceSum(amaxSf_)().primitiveField()),
              CoNum_(0),
              meanCoNum_(0)
        {
            // Only need to evaluate this once at construction time
            // Avoid evaluation every time update is called
            if (mesh_.nInternalFaces() == false)
            {
                throw std::invalid_argument("No internal mesh faces.");
            }
        };

        // Destructor
        inline ~acousticCourantNo(){};

        // Calculate and print acoustic Courant number
        inline void calculate()
        {
            amaxSf_ = 0.5 * (mag((fvc::interpolate(U_) & mesh_.Sf())) + (fvc::interpolate(sqrt(thermo_.gamma() * (thermo_.Cp() - thermo_.Cv()) * thermo_.T())) * mesh_.magSf()));
            sumAmaxSf_ = fvc::surfaceSum(amaxSf_)().primitiveField();
            CoNum_ = 0.5 * gMax(sumAmaxSf_ / mesh_.V().field()) * runTime_.deltaTValue();
            meanCoNum_ = 0.5 * (gSum(sumAmaxSf_) / gSum(mesh_.V().field())) * runTime_.deltaTValue();

            // Output
            Info << "Mean and max acoustic Courant Numbers = " << meanCoNum_ << " " << CoNum_ << endl;
        };

        // Update the time step
        inline void setDeltaT(
            const scalar maxDeltaT,
            const scalar maxCo,
            const bool adjustTimeStep)
        {
            if (adjustTimeStep)
            {
                const scalar maxDeltaTFact = maxCo / (CoNum_ + SMALL);
                const scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2);

                runTime_.setDeltaT(
                    min(
                        deltaTFact * runTime_.deltaTValue(),
                        maxDeltaT));

                Info << "deltaT = " << runTime_.deltaTValue() << endl;
            }
        }

    private:
        // Const access to velocity field
        const volVectorField &U_;

        // Const access to mesh
        const fvMesh &mesh_;

        // Const access to thermo
        const Thermo &thermo_;

        // Non-const access to runTime
        Time &runTime_;

        surfaceScalarField amaxSf_;
        scalarField sumAmaxSf_;

        // Max Courant number
        scalar CoNum_;

        // Mean Courant number
        scalar meanCoNum_;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
