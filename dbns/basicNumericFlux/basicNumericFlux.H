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
    Foam::basicNumericFlux

Description
    Base class for runtime-selectable flux and limiters

Author
    Henrik Rusche

SourceFiles
    basicNumericFlux.C
    newBasicNumericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef basicNumericFlux_H
#define basicNumericFlux_H

#include "fvMesh.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                             Class basicNumericFlux Declaration
    \*---------------------------------------------------------------------------*/

    class basicNumericFlux
    {
        // Private data

        //- Reference to mesh
        const fvMesh &mesh_;

    protected:
        // Protected member functions

        //- Construct as copy (not implemented)
        basicNumericFlux(const basicNumericFlux &);

    public:
        //- Runtime type information
        TypeName("basicNumericFlux");

        //- Declare run-time constructor selection table

#ifndef SWIG
        declareRunTimeSelectionTable(
            autoPtr,
            basicNumericFlux,
            state,
            (
                const volScalarField &p,
                const volVectorField &U,
                const volScalarField &T,
                const volScalarField &upwindingFactor,
                const basicThermo &thermo),
            (p, U, T, upwindingFactor, thermo));
#endif

        // Constructors

        //- Construct from components
        basicNumericFlux(
            const fvMesh &mesh);

        //- Selector
        static autoPtr<basicNumericFlux> New(
            const volScalarField &p,
            const volVectorField &U,
            const volScalarField &T,
            const volScalarField &upwindingFactor,
            const basicThermo &thermo);

        //- Destructor
        virtual ~basicNumericFlux();

        // Member functions

        //- Return mesh reference
        const fvMesh &mesh() const
        {
            return mesh_;
        }

        // Return fluxes

        //- Return density flux
        virtual const surfaceScalarField &rhoFlux() const = 0;

        //- Return velocity flux
        virtual const surfaceVectorField &rhoUFlux() const = 0;

        //- Return energy flux
        virtual const surfaceScalarField &rhoEFlux() const = 0;

        // Update fluxes based on current state

        //- Compute flux
        virtual void computeFlux() = 0;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
