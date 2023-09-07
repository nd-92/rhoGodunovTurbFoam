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
    Density-based shock-capturing compressible flow solver based on the Roe
    flux difference splitting scheme of Bui with absorbing sponge zones
    and explicitly bounded fields

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

// #include "dimensionedScalar.H"
// #include "volFieldsFwd.H"
// #include "volFields.H"
// #include "fvc.H"

// void boundMinMax(
//     volScalarField &vsf,
//     const dimensionedScalar &lowerBound,
//     const dimensionedScalar &upperBound)
// {
//     const scalar minVsf = min(vsf).value();
//     const scalar maxVsf = max(vsf).value();

//     if (minVsf < lowerBound.value() || maxVsf > upperBound.value())
//     {
//         Info << "bounding " << vsf.name()
//              << ", min: " << minVsf
//              << " max: " << maxVsf
//              << " average: " << gAverage(vsf.internalField())
//              << endl;
//     }

//     // if (minVsf < lowerBound.value())
//     // {
//     //     vsf.primitiveFieldRef() = max(
//     //         max(
//     //             vsf.primitiveField(),
//     //             fvc::average(max(vsf, lowerBound))().primitiveField() * pos0(lowerBound.value() - vsf.primitiveField())),
//     //         lowerBound.value());
//     //     Info << "new max: " << gMax(vsf.internalField()) << endl;
//     //     vsf.boundaryFieldRef() = max(vsf.boundaryField(), lowerBound.value());
//     // }

//     // if (maxVsf > upperBound.value())
//     // {
//     //     vsf.primitiveFieldRef() = min(
//     //         min(
//     //             vsf.primitiveField(),
//     //             fvc::average(min(vsf, upperBound))().primitiveField() * neg(upperBound.value() - vsf.primitiveField())),
//     //         upperBound.value());
//     //     Info << "new max: " << gMax(vsf.internalField()) << endl;
//     //     vsf.boundaryFieldRef() = min(vsf.boundaryField(), lowerBound.value());
//     // }

//     if (minVsf < lowerBound.value())
//     {
//         vsf.primitiveFieldRef() = max(
//             max(
//                 vsf.primitiveField(),
//                 fvc::average(max(vsf, lowerBound))().primitiveField() * pos0(lowerBound.value() - vsf.primitiveField())),
//             lowerBound.value());
//         Info << "new max: " << gMax(vsf.internalField()) << endl;
//         vsf.boundaryFieldRef() = max(vsf.boundaryField(), lowerBound.value());
//     }

//     if (maxVsf > upperBound.value())
//     {
//         vsf.primitiveFieldRef() = min(
//             min(
//                 vsf.primitiveField(),
//                 fvc::average(min(vsf, upperBound))().primitiveField() * neg(upperBound.value() - vsf.primitiveField())),
//             upperBound.value());
//         Info << "new max: " << gMax(vsf.internalField()) << endl;
//         vsf.boundaryFieldRef() = min(vsf.boundaryField(), lowerBound.value());
//     }

//     // return vsf;
// }

int main(int argc, char *argv[])
{

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
#include "createTimeControls.H"
    // #include "readFieldBounds.H"

    Info << "Starting time loop" << endl;

    if (applyDamping == true)
    {
#include "dampingFields.H"
        while (runTime.run())
        {
            // Execute main solver loop
#include "buiTurbFoam.H"

            // Apply field bounds
            // boundMinMax(he, heMin, heMax);
            // boundMinMax(rho, rhoMin, rhoMax);
            // boundMinMax(T, TMin, TMax);
            // boundMinMax(p, pMin, pMax);

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

        Info << "\nEnd\n"
             << endl;

        return 0;
    }
    else
    {
        while (runTime.run())
        {
            // Execute main solver loop
#include "buiTurbFoam.H"

            // Correct turbulence fields
            turbulence->correct();

            // Write runtime output
            runTime.write();
            runTime.printExecutionTime(Info);
        }

        Info << "\nEnd\n"
             << endl;

        return 0;
    }
}

// ************************************************************************* //
