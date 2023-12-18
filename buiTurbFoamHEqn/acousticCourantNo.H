/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
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

Global
    acousticCourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

surfaceScalarField amaxSf(                                                                                       //
    0.5 *                                                                                                        //
    (mag(((fvc::interpolate(rhoU) / fvc::interpolate(rho)) & mesh.Sf()))                                         //
     + (fvc::interpolate(sqrt((thermo.Cp() / thermo.Cv()) * (thermo.Cp() - thermo.Cv()) * T)) * mesh.magSf()))); //

if (mesh.nInternalFaces())
{
    scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    CoNum = 0.5 * gMax(sumAmaxSf / mesh.V().field()) * runTime.deltaTValue();

    meanCoNum = 0.5 * (gSum(sumAmaxSf) / gSum(mesh.V().field())) * runTime.deltaTValue();
}

Info << "Mean and max acoustic Courant Numbers = " << meanCoNum << " " << CoNum << endl;

// ************************************************************************* //
