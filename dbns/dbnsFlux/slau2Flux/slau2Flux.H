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
    slau2Flux

Description
    Implementation of SLAU2 approximate Riemann solver scheme

SourceFiles
    buiFlux.H

\*---------------------------------------------------------------------------*/

#ifndef slau2Flux_H
#define slau2Flux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                               Class slau2Flux Declaration
    \*---------------------------------------------------------------------------*/

    class slau2Flux
    {
    public:
        // Constructors

        //- Construct null
        slau2Flux()
        {
        }

        // Destructor - default

        // Member Functions

        //- Update fluxes based on current left and right states
        void evaluateFlux(
            scalar &rhoFlux,
            vector &rhoUFlux,
            scalar &rhoEFlux,
            const scalar &pLeft,
            const scalar &pRight,
            const vector &ULeft,
            const vector &URight,
            const scalar &TLeft,
            const scalar &TRight,
            const scalar &RLeft,
            const scalar &RRight,
            const scalar &CvLeft,
            const scalar &CvRight,
            const vector &Sf,
            const scalar &magSf,
            const scalar &buiEps) const;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
