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
    rhoCentralFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver based on
    central-upwind schemes of Kurganov and Tadmor with
    support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "maxwellSlipUFvPatchVectorField.H"
#include "smoluchowskiJumpTFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote(
        "Density-based compressible flow solver based on"
        " central-upwind schemes of Kurganov and Tadmor with"
        " support for mesh-motion and topology changes.");

#define NO_CONTROL
#include "postProcess.H"
#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
// #include "createDynamicFvMesh.H"
#include "createFields.H"
#include "createFieldRefs.H"
#include "createTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const dimensionedScalar v_zero(dimVolume / dimTime, Zero);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info << endl;
    Info << "Starting time loop" << endl;
    Info << endl;

    while (runTime.run())
    {

        // Execute main body of solver
#include "readTimeControls.H"

        if (!LTS)
        {
#include "setDeltaT.H"

            ++runTime;

            // Do any mesh changes
            // mesh.update();
        }

        // --- Directed interpolation of primitive fields onto faces

        const surfaceScalarField rho_pos(interpolate(rho, pos));
        const surfaceScalarField rho_neg(interpolate(rho, neg));

        const surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        const surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        const volScalarField rPsi("rPsi", 1.0 / psi);
        const surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        const surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        const surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        const surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        const surfaceVectorField U_pos("U_pos", rhoU_pos / rho_pos);
        const surfaceVectorField U_neg("U_neg", rhoU_neg / rho_neg);

        const surfaceScalarField p_pos("p_pos", rho_pos * rPsi_pos);
        const surfaceScalarField p_neg("p_neg", rho_neg * rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        // Note: extracted out the orientation so becomes unoriented
        phiv_pos.setOriented(false);
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
        phiv_neg.setOriented(false);

        // Make fluxes relative to mesh-motion
        // if (mesh.moving())
        // {
        //     surfaceScalarField meshPhi(mesh.phi());
        //     meshPhi.setOriented(false);
        //     phiv_pos -= meshPhi;
        //     phiv_neg -= meshPhi;
        // }

        const volScalarField c("c", sqrt(thermo.Cp() / thermo.Cv() * rPsi));

        const surfaceScalarField cSf_pos("cSf_pos", interpolate(c, pos, T.name()) * mesh.magSf());

        const surfaceScalarField cSf_neg("cSf_neg", interpolate(c, neg, T.name()) * mesh.magSf());

        const surfaceScalarField ap("ap", max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero));

        const surfaceScalarField am("am", min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero));

        const surfaceScalarField a_pos("a_pos", ap / (ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        const surfaceScalarField aSf("aSf", am * a_pos);

        const surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        const surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        const surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        if (mesh.nInternalFaces())
        {
            const scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

            CoNum = 0.5 * gMax(sumAmaxSf / mesh.V().field()) * runTime.deltaTValue();

            meanCoNum = 0.5 * (gSum(sumAmaxSf) / gSum(mesh.V().field())) * runTime.deltaTValue();

            Info << "Mean and max Courant Numbers = " << meanCoNum << " " << CoNum << endl;
        }

        if (LTS)
        {
            {
                volScalarField &rDeltaT = trDeltaT.ref();

                const scalar rDeltaTSmoothingCoeff(runTime.controlDict().getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02));

                // Set the reciprocal time-step from the local Courant number
                rDeltaT.ref() = max(1 / dimensionedScalar("maxDeltaT", dimTime, maxDeltaT), fvc::surfaceSum(amaxSf)()() / ((2 * maxCo) * mesh.V()));

                // Update tho boundary values of the reciprocal time-step
                rDeltaT.correctBoundaryConditions();

                fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

                Info << "Flow time scale min/max = " << gMin(1 / rDeltaT.primitiveField()) << ", " << gMax(1 / rDeltaT.primitiveField()) << endl;
            }

            ++runTime;
        }

        Info << "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos * rho_pos + aphiv_neg * rho_neg;

        surfaceVectorField phiU(aphiv_pos * rhoU_pos + aphiv_neg * rhoU_neg);
        // Note: reassembled orientation from the pos and neg parts so becomes
        // oriented
        phiU.setOriented(true);

        const surfaceVectorField phiUp(phiU + (a_pos * p_pos + a_neg * p_neg) * mesh.Sf());

        const surfaceScalarField phiEp("phiEp", aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg);

        // Make flux for pressure-work absolute
        // if (mesh.moving())
        // {
        //     surfaceScalarField meshPhi(mesh.phi());
        //     meshPhi.setOriented(false);
        //     phiEp += meshPhi * (a_pos * p_pos + a_neg * p_neg);
        // }

        const volScalarField muEff("muEff", turbulence->muEff());
        const volTensorField tauMC("tauMC", muEff * dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() = rhoU() / rho();
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

        if (!inviscid)
        {
            solve(                         //
                fvm::ddt(rho, U)           //
                - fvc::ddt(rho, U)         //
                - fvm::laplacian(muEff, U) //
                - fvc::div(tauMC));        //
            rhoU = rho * U;                //
        }

        // --- Solve energy
        const surfaceScalarField sigmaDotU(                          //
            "sigmaDotU",                                             //
            (fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) //
             + fvc::dotInterpolate(mesh.Sf(), tauMC)) &              //
                (a_pos * U_pos + a_neg * U_neg));                    //

        solve(                      //
            fvm::ddt(rhoE)          //
            + fvc::div(phiEp)       //
            - fvc::div(sigmaDotU)); //

        e = rhoE / rho - 0.5 * magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

        if (!inviscid)
        {
            solve(                                            //
                fvm::ddt(rho, e)                              //
                - fvc::ddt(rho, e)                            //
                - fvm::laplacian(turbulence->alphaEff(), e)); //
            thermo.correct();                                 //
            rhoE = rho * (e + 0.5 * magSqr(U));               //
        }

        p.ref() = rho() / psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

        // Apply damping
        T = (T_inf * acousticDamping) + (T * noDamping);
        T.correctBoundaryConditions();
        p = (p_inf * acousticDamping) + (p * noDamping);
        p.correctBoundaryConditions();
        rho = (rho_inf * acousticDamping) + (rho * noDamping);
        rho.correctBoundaryConditions();
        U = (U_inf * acousticDamping) + (U * noDamping);
        U.correctBoundaryConditions();

        turbulence->correct();

        // Do I/O
        runTime.write();

        // Print execution time
        runTime.printExecutionTime(Info);
    }

    Info << "End" << endl;

    return 0;
}

// ************************************************************************* //
