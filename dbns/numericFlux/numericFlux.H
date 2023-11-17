/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    numericFlux

Description
    Single level numeric flux class for density-based solvers

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    numericFlux.H
    numericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef numericFlux_H
#define numericFlux_H

#include "numericFluxBase.H"
#include "basicThermo.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                              Class numericFlux Declaration
    \*---------------------------------------------------------------------------*/

    template <class Flux, class Limiter>
    class numericFlux
        : public numericFluxBase<Flux>
    {
        // Private data

        // Reference to primitive fields

        //- Static pressure
        const volScalarField &p_;

        //- Velocity
        const volVectorField &U_;

        //- Static temperature
        const volScalarField &T_;

        //- Upwinding factor
        const surfaceScalarField buiEps_;

        //- Reference to the thermophysicalModel
        const basicThermo &thermo_;

        // Fluxes

        //- Density flux
        surfaceScalarField rhoFlux_;

        //- Velocity flux
        surfaceVectorField rhoUFlux_;

        //- Energy flux
        surfaceScalarField rhoEFlux_;

        // Gradients

        //- Static pressure gradient
        volVectorField gradP;

        //- Velocity gradient
        volTensorField gradU;

        //- Static temperature gradient
        volVectorField gradT;

        // Private Member Functions

        //- Disallow default bitwise copy construct
        numericFlux(const numericFlux &);

        //- Disallow default bitwise assignment
        void operator=(const numericFlux &);

        //- Return internal field of mass flux
        const scalarField &rhoFluxI() const
        {
            return rhoFlux_.internalField();
        }

        //- Return internal field of momentum flux
        const vectorField &rhoUFluxI() const
        {
            return rhoUFlux_.internalField();
        }

        //- Return access to internal field of energy flux
        const scalarField &rhoEFluxI() const
        {
            return rhoEFlux_.internalField();
        }

    public:
        //- Runtime type information
        TypeName("numericFlux");

        // Constructors

        //- Construct from components
        numericFlux(
            const volScalarField &p,
            const volVectorField &U,
            const volScalarField &T,
            const volScalarField &upwindingFactor,
            const basicThermo &thermo);

        //- Destructor
        virtual ~numericFlux()
        {
        }

        // Member Functions

        // Return fluxes

        //- Return density flux
        virtual const surfaceScalarField &rhoFlux() const
        {
            return rhoFlux_;
        }

        //- Return velocity flux
        virtual const surfaceVectorField &rhoUFlux() const
        {
            return rhoUFlux_;
        }

        //- Return energy flux
        virtual const surfaceScalarField &rhoEFlux() const
        {
            return rhoEFlux_;
        }

        // Update fluxes based on current state

        //- Compute flux
        virtual void computeFlux();
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "numericFlux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
