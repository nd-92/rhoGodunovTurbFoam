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
#include "MDLimiter.H"
#include "cyclicFvPatch.H"
#include "wallFvPatch.H"
#include "shockSensor.H"
#include "fvCFD.H"
#include "buiFlux.H"
#include "slau2Flux.H"
#include "slauFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                              Class numericFlux Declaration
    \*---------------------------------------------------------------------------*/

    template <class Flux, class Limiter>
    class numericFlux : public numericFluxBase<Flux>
    {
        // Private data

        // Reference to primitive fields

        //- Static pressure
        const volScalarField &p_;

        //- Velocity
        const volVectorField &U_;

        //- Static temperature
        const volScalarField &T_;

        //- Upwinding factor vol vield
        const volScalarField &upwindingFactor_;

        //- Upwinding factor surface field
        surfaceScalarField buiEps_;

        //- Minimum upwinding factor
        const scalar upwindingFactorMin_;

        //- Reference to the thermophysicalModel
        const basicThermo &thermo_;

        //-Thermodynamic properties
        const volScalarField Cv_;
        volScalarField R_;

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
        inline numericFlux(const numericFlux &);

        //- Disallow default bitwise assignment
        inline void operator=(const numericFlux &);

        //- Return internal field of mass flux
        [[nodiscard]] inline const scalarField &rhoFluxI() const
        {
            return rhoFlux_.internalField();
        }

        //- Return internal field of momentum flux
        [[nodiscard]] inline const vectorField &rhoUFluxI() const
        {
            return rhoUFlux_.internalField();
        }

        //- Return access to internal field of energy flux
        [[nodiscard]] inline const scalarField &rhoEFluxI() const
        {
            return rhoEFlux_.internalField();
        }

    public:
        //- Runtime type information
        inline TypeName("numericFlux");

        // Constructors

        //- Construct from components
        inline numericFlux(
            const volScalarField &p,
            const volVectorField &U,
            const volScalarField &T,
            const volScalarField &upwindingFactor,
            const basicThermo &thermo)
            : numericFluxBase<Flux>(p.mesh()),
              p_(p),
              U_(U),
              T_(T),
              upwindingFactor_(upwindingFactor),
              buiEps_(fvc::interpolate(upwindingFactor)),
              upwindingFactorMin_(this->mesh().schemesDict().subDict("divSchemes").subDict("dbns").lookupOrDefault("epsilonMin", 1.0)),
              thermo_(thermo),
              Cv_(thermo.Cv()),
              R_(thermo.Cp() - thermo.Cv()),
              rhoFlux_(
                  IOobject(
                      "phi",
                      this->mesh().time().timeName(),
                      this->mesh(),
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                  (fvc::interpolate(thermo.rho() * U) & this->mesh().Sf())),
              rhoUFlux_(
                  IOobject(
                      "rhoUFlux",
                      this->mesh().time().timeName(),
                      this->mesh(),
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                  rhoFlux_ * fvc::interpolate(U)),
              rhoEFlux_(
                  IOobject(
                      "rhoEFlux",
                      this->mesh().time().timeName(),
                      this->mesh(),
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                  rhoFlux_ * fvc::interpolate(thermo.Cv() * T + 0.5 * magSqr(U))),
              gradP(fvc::grad(p)),
              gradU(fvc::grad(U)),
              gradT(fvc::grad(T)){};

        //- Destructor
        inline virtual ~numericFlux()
        {
        }

        // Member Functions

        // Return fluxes

        //- Return density flux
        [[nodiscard]] virtual inline const surfaceScalarField &rhoFlux() const
        {
            return rhoFlux_;
        }

        //- Return velocity flux
        [[nodiscard]] virtual inline const surfaceVectorField &rhoUFlux() const
        {
            return rhoUFlux_;
        }

        //- Return energy flux
        [[nodiscard]] virtual inline const surfaceScalarField &rhoEFlux() const
        {
            return rhoEFlux_;
        }

        // Update fluxes based on current state

        //- Compute flux
        // inline void computeFlux()
        // {
        //     // Update the upwinding factor
        //     buiEps_ = (linearInterpolate(upwindingFactor_) * (1 - upwindingFactorMin_)) + upwindingFactorMin_;

        //     // Get face-to-cell addressing: face area point from owner to neighbour
        //     const labelUList &owner = this->mesh().owner();
        //     const labelUList &neighbour = this->mesh().neighbour();

        //     // Get the face area vector
        //     const surfaceVectorField &Sf = this->mesh().Sf();
        //     const surfaceScalarField &magSf = this->mesh().magSf();

        //     const volVectorField &cellCentre = this->mesh().C();
        //     const surfaceVectorField &faceCentre = this->mesh().Cf();

        //     // Get gradients
        //     // Coupled patch update on gradients moved into gradScheme.C
        //     // HJ, 22/Apr/2016;

        //     // Changed return type for gradient cacheing.  HJ, 22/Apr/2016
        //     const tmp<volVectorField> tgradP = fvc::grad(p_);
        //     const volVectorField &gradP = tgradP();

        //     const tmp<volTensorField> tgradU = fvc::grad(U_);
        //     const volTensorField &gradU = tgradU();

        //     const tmp<volVectorField> tgradT = fvc::grad(T_);
        //     const volVectorField &gradT = tgradT();

        //     const MDLimiter<scalar, Limiter> scalarPLimiter(
        //         this->p_,
        //         gradP);

        //     const MDLimiter<vector, Limiter> vectorULimiter(
        //         this->U_,
        //         gradU);

        //     const MDLimiter<scalar, Limiter> scalarTLimiter(
        //         this->T_,
        //         gradT);

        //     // Get limiters
        //     const volScalarField &pLimiter = scalarPLimiter.phiLimiter();
        //     const volVectorField &ULimiter = vectorULimiter.phiLimiter();
        //     const volScalarField &TLimiter = scalarTLimiter.phiLimiter();

        //     // Calculate fluxes at internal faces
        //     forAll(owner, faceI)
        //     {
        //         const label own = owner[faceI];
        //         const label nei = neighbour[faceI];

        //         const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
        //         const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

        //         // calculate fluxes with reconstructed primitive variables at faces
        //         Flux::evaluateFlux(
        //             rhoFlux_[faceI],
        //             rhoUFlux_[faceI],
        //             rhoEFlux_[faceI],
        //             p_[own] + pLimiter[own] * (deltaRLeft & gradP[own]),
        //             p_[nei] + pLimiter[nei] * (deltaRRight & gradP[nei]),
        //             U_[own] + cmptMultiply(ULimiter[own], (deltaRLeft & gradU[own])),
        //             U_[nei] + cmptMultiply(ULimiter[nei], (deltaRRight & gradU[nei])),
        //             T_[own] + TLimiter[own] * (deltaRLeft & gradT[own]),
        //             T_[nei] + TLimiter[nei] * (deltaRRight & gradT[nei]),
        //             R_[own],
        //             R_[nei],
        //             Cv_[own],
        //             Cv_[nei],
        //             Sf[faceI],
        //             magSf[faceI],
        //             buiEps_[faceI]);
        //         // Flux::evaluateFlux(
        //         //     rhoFlux_[faceI],
        //         //     rhoUFlux_[faceI],
        //         //     rhoEFlux_[faceI],
        //         //     p_[own] + (deltaRLeft & gradP[own]),
        //         //     p_[nei] + (deltaRRight & gradP[nei]),
        //         //     U_[own] + (deltaRLeft & gradU[own]),
        //         //     U_[nei] + (deltaRRight & gradU[nei]),
        //         //     T_[own] + (deltaRLeft & gradT[own]),
        //         //     T_[nei] + (deltaRRight & gradT[nei]),
        //         //     R_[own],
        //         //     R_[nei],
        //         //     Cv_[own],
        //         //     Cv_[nei],
        //         //     Sf[faceI],
        //         //     magSf[faceI],
        //         //     // buiEps_[faceI]);
        //         //     1.0);
        //     }

        //     // Update boundary field and values
        //     forAll(rhoFlux_.boundaryField(), patchi)
        //     {
        //         const fvPatch &curPatch = p_.boundaryField()[patchi].patch();

        //         // Fluxes
        //         fvsPatchScalarField &pRhoFlux = rhoFlux_.boundaryFieldRef()[patchi];
        //         fvsPatchVectorField &pRhoUFlux = rhoUFlux_.boundaryFieldRef()[patchi];
        //         fvsPatchScalarField &pRhoEFlux = rhoEFlux_.boundaryFieldRef()[patchi];

        //         // Patch fields
        //         const fvPatchScalarField &pp = p_.boundaryField()[patchi];
        //         const vectorField &pU = U_.boundaryField()[patchi];
        //         const scalarField &pT = T_.boundaryField()[patchi];

        //         const scalarField &pCv = Cv_.boundaryField()[patchi];
        //         const scalarField &pR = R_.boundaryField()[patchi];

        //         // Gradients
        //         const fvPatchVectorField &pGradP = gradP.boundaryField()[patchi];
        //         const fvPatchTensorField &pGradU = gradU.boundaryField()[patchi];
        //         const fvPatchVectorField &pGradT = gradT.boundaryField()[patchi];

        //         // Limiters
        //         const fvPatchScalarField &pPatchLim = pLimiter.boundaryField()[patchi];
        //         const fvPatchVectorField &UPatchLim = ULimiter.boundaryField()[patchi];
        //         const fvPatchScalarField &TPatchLim = TLimiter.boundaryField()[patchi];

        //         // Face areas
        //         const fvsPatchVectorField &pSf = Sf.boundaryField()[patchi];
        //         const fvsPatchScalarField &pMagSf = magSf.boundaryField()[patchi];

        //         if (pp.coupled())
        //         {
        //             // Coupled patch
        //             const scalarField ppLeft = p_.boundaryField()[patchi].patchInternalField();
        //             const scalarField ppRight = p_.boundaryField()[patchi].patchNeighbourField();

        //             const vectorField pULeft = U_.boundaryField()[patchi].patchInternalField();
        //             const vectorField pURight = U_.boundaryField()[patchi].patchNeighbourField();

        //             const scalarField pTLeft = T_.boundaryField()[patchi].patchInternalField();
        //             const scalarField pTRight = T_.boundaryField()[patchi].patchNeighbourField();

        //             // Gradients
        //             const vectorField pgradPLeft = pGradP.patchInternalField();
        //             const vectorField pgradPRight = pGradP.patchNeighbourField();

        //             const tensorField pgradULeft = pGradU.patchInternalField();
        //             const tensorField pgradURight = pGradU.patchNeighbourField();

        //             const vectorField pgradTLeft = pGradT.patchInternalField();
        //             const vectorField pgradTRight = pGradT.patchNeighbourField();

        //             // Geometry: call the raw cell-to-face vector by calling
        //             // the base patch (cell-to-face) delta coefficient
        //             // Work out the right delta from the cell-to-cell delta
        //             // across the coupled patch and left delta
        //             const vectorField pDeltaRLeft = curPatch.fvPatch::delta();
        //             const vectorField pDeltaRRight = pDeltaRLeft - curPatch.delta();

        //             // Limiters

        //             const scalarField ppLimiterLeft = pPatchLim.patchInternalField();
        //             const scalarField ppLimiterRight = pPatchLim.patchNeighbourField();

        //             const vectorField pULimiterLeft = UPatchLim.patchInternalField();
        //             const vectorField pULimiterRight = UPatchLim.patchNeighbourField();

        //             const scalarField pTLimiterLeft = TPatchLim.patchInternalField();
        //             const scalarField pTLimiterRight = TPatchLim.patchNeighbourField();

        //             forAll(pp, facei)
        //             {
        //                 Flux::evaluateFlux(
        //                     pRhoFlux[facei],
        //                     pRhoUFlux[facei],
        //                     pRhoEFlux[facei],
        //                     ppLeft[facei] + ppLimiterLeft[facei] * (pDeltaRLeft[facei] & pgradPLeft[facei]),
        //                     ppRight[facei] + ppLimiterRight[facei] * (pDeltaRRight[facei] & pgradPRight[facei]),
        //                     pULeft[facei] + cmptMultiply(pULimiterLeft[facei], pDeltaRLeft[facei] & pgradULeft[facei]),
        //                     pURight[facei] + cmptMultiply(pULimiterRight[facei], pDeltaRRight[facei] & pgradURight[facei]),
        //                     pTLeft[facei] + pTLimiterLeft[facei] * (pDeltaRLeft[facei] & pgradTLeft[facei]),
        //                     pTRight[facei] + pTLimiterRight[facei] * (pDeltaRRight[facei] & pgradTRight[facei]),
        //                     pR[facei],
        //                     pR[facei],
        //                     pCv[facei],
        //                     pCv[facei],
        //                     pSf[facei],
        //                     pMagSf[facei],
        //                     buiEps_[facei]);
        //                 // Flux::evaluateFlux(
        //                 //     pRhoFlux[facei],
        //                 //     pRhoUFlux[facei],
        //                 //     pRhoEFlux[facei],
        //                 //     ppLeft[facei] + (pDeltaRLeft[facei] & pgradPLeft[facei]),
        //                 //     ppRight[facei] + (pDeltaRRight[facei] & pgradPRight[facei]),
        //                 //     pULeft[facei] + pDeltaRLeft[facei] & pgradULeft[facei],
        //                 //     pURight[facei] + pDeltaRRight[facei] & pgradURight[facei],
        //                 //     pTLeft[facei] + (pDeltaRLeft[facei] & pgradTLeft[facei]),
        //                 //     pTRight[facei] + (pDeltaRRight[facei] & pgradTRight[facei]),
        //                 //     pR[facei],
        //                 //     pR[facei],
        //                 //     pCv[facei],
        //                 //     pCv[facei],
        //                 //     pSf[facei],
        //                 //     pMagSf[facei],
        //                 //     // buiEps_[facei]);
        //                 //     1.0);
        //             }
        //         }
        //         else
        //         {
        //             // Coupled patch
        //             // const scalarField pp = p_.boundaryField()[patchi].patchInternalField();
        //             // const scalarField ppRight = p_.boundaryField()[patchi].patchNeighbourField();

        //             // const vectorField pU = U_.boundaryField()[patchi].patchInternalField();
        //             // const vectorField pURight = U_.boundaryField()[patchi].patchNeighbourField();

        //             // const scalarField pT = T_.boundaryField()[patchi].patchInternalField();
        //             // const scalarField pT = T_.boundaryField()[patchi].patchNeighbourField();

        //             // Gradients
        //             const vectorField pgradP = pGradP.patchInternalField();
        //             // const vectorField pgradPRight = pGradP.patchNeighbourField();

        //             const tensorField pgradU = pGradU.patchInternalField();
        //             // const tensorField pgradURight = pGradU.patchNeighbourField();

        //             const vectorField pgradT = pGradT.patchInternalField();
        //             // const vectorField pgradTRight = pGradT.patchNeighbourField();

        //             // Geometry: call the raw cell-to-face vector by calling
        //             // the base patch (cell-to-face) delta coefficient
        //             // Work out the right delta from the cell-to-cell delta
        //             // across the coupled patch and left delta
        //             const vectorField pDeltaRLeft = curPatch.fvPatch::delta();
        //             const vectorField pDeltaRRight = pDeltaRLeft - curPatch.delta();

        //             // Limiters

        //             const scalarField ppLimiter = pPatchLim.patchInternalField();
        //             // const scalarField ppLimiterRight = pPatchLim.patchNeighbourField();

        //             const vectorField pULimiter = UPatchLim.patchInternalField();
        //             // const vectorField pULimiterRight = UPatchLim.patchNeighbourField();

        //             const scalarField pTLimiter = TPatchLim.patchInternalField();
        //             // const scalarField pTLimiterRight = TPatchLim.patchNeighbourField();

        //             forAll(pp, facei)
        //             {
        //                 // Calculate fluxes
        //                 Flux::evaluateFlux(
        //                     pRhoFlux[facei],
        //                     pRhoUFlux[facei],
        //                     pRhoEFlux[facei],

        //                     pp[facei] + ppLimiter[facei] * (pDeltaRLeft[facei] & pgradP[facei]),
        //                     pp[facei] + ppLimiter[facei] * (pDeltaRRight[facei] & pgradP[facei]),
        //                     pU[facei] + cmptMultiply(pULimiter[facei], pDeltaRLeft[facei] & pgradU[facei]),
        //                     pU[facei] + cmptMultiply(pULimiter[facei], pDeltaRRight[facei] & pgradU[facei]),
        //                     pT[facei] + pTLimiter[facei] * (pDeltaRLeft[facei] & pgradT[facei]),
        //                     pT[facei] + pTLimiter[facei] * (pDeltaRRight[facei] & pgradT[facei]),

        //                     // pp[facei],
        //                     // pp[facei],
        //                     // pU[facei],
        //                     // pU[facei],
        //                     // pT[facei],
        //                     // pT[facei],

        //                     pR[facei],
        //                     pR[facei],
        //                     pCv[facei],
        //                     pCv[facei],
        //                     pSf[facei],
        //                     pMagSf[facei],
        //                     buiEps_[facei]);
        //                 // 1.0);
        //             }
        //         }
        //     }
        // }

        // template <class Flux, class Limiter>
        inline void computeFlux()
        {
            // Get face-to-cell addressing: face area point from owner to neighbour
            const labelUList &owner = this->mesh().owner();
            const labelUList &neighbour = this->mesh().neighbour();

            // Get the face area vector
            const surfaceVectorField &Sf = this->mesh().Sf();
            const surfaceScalarField &magSf = this->mesh().magSf();

            const volVectorField &cellCentre = this->mesh().C();
            const surfaceVectorField &faceCentre = this->mesh().Cf();

            // Thermodynamics
            const volScalarField Cv = thermo_.Cv();
            const volScalarField R = thermo_.Cp() - Cv;

            gradP = fvc::grad(p_);
            gradP.correctBoundaryConditions();

            gradU = fvc::grad(U_);
            gradU.correctBoundaryConditions();

            gradT = fvc::grad(T_);
            gradT.correctBoundaryConditions();

            // Calculate fluxes at internal faces
            forAll(owner, faceI)
            {
                const label own = owner[faceI];
                const label nei = neighbour[faceI];

                const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
                const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

                // calculate fluxes with reconstructed primitive variables at faces
                Flux::evaluateFlux(
                    rhoFlux_[faceI],
                    rhoUFlux_[faceI],
                    rhoEFlux_[faceI],
                    p_[own] + (deltaRLeft & gradP[own]),
                    p_[nei] + (deltaRRight & gradP[nei]),
                    U_[own] + (deltaRLeft & gradU[own]),
                    U_[nei] + (deltaRRight & gradU[nei]),
                    T_[own] + (deltaRLeft & gradT[own]),
                    T_[nei] + (deltaRRight & gradT[nei]),
                    R[own],
                    R[nei],
                    Cv[own],
                    Cv[nei],
                    Sf[faceI],
                    magSf[faceI],
                    buiEps_[faceI]);
            }

            // Update boundary field and values
            forAll(p_.boundaryField(), patchi)
            {
                const fvPatch &curPatch = p_.boundaryField()[patchi].patch();

                // Fluxes
                fvsPatchScalarField &pRhoFlux = rhoFlux_.boundaryFieldRef()[patchi];
                fvsPatchVectorField &pRhoUFlux = rhoUFlux_.boundaryFieldRef()[patchi];
                fvsPatchScalarField &pRhoEFlux = rhoEFlux_.boundaryFieldRef()[patchi];

                // Patch fields
                const fvPatchScalarField &pp = p_.boundaryField()[patchi];
                const vectorField &pU = U_.boundaryField()[patchi];
                const scalarField &pT = T_.boundaryField()[patchi];

                const scalarField &pCv = Cv.boundaryField()[patchi];
                const scalarField &pR = R.boundaryField()[patchi];

                // Gradients
                const fvPatchVectorField &pGradP = gradP.boundaryField()[patchi];
                const fvPatchTensorField &pGradU = gradU.boundaryField()[patchi];
                const fvPatchVectorField &pGradT = gradT.boundaryField()[patchi];

                // Face areas
                const fvsPatchVectorField &pSf = Sf.boundaryField()[patchi];
                const fvsPatchScalarField &pMagSf = magSf.boundaryField()[patchi];
                const fvsPatchScalarField &pBuiEps = buiEps_.boundaryField()[patchi];

                const fvPatchVectorField &pCellCenter = cellCentre.boundaryField()[patchi];

                if (pp.coupled())
                {
                    // Coupled patch
                    const scalarField ppLeft =
                        p_.boundaryField()[patchi].patchInternalField();

                    const scalarField ppRight =
                        p_.boundaryField()[patchi].patchNeighbourField();

                    const vectorField pULeft =
                        U_.boundaryField()[patchi].patchInternalField();

                    const vectorField pURight =
                        U_.boundaryField()[patchi].patchNeighbourField();

                    const scalarField pTLeft =
                        T_.boundaryField()[patchi].patchInternalField();

                    const scalarField pTRight =
                        T_.boundaryField()[patchi].patchNeighbourField();

                    // Gradients
                    const vectorField pgradPLeft = pGradP.patchInternalField();
                    const vectorField pgradPRight = pGradP.patchNeighbourField();

                    const tensorField pgradULeft = pGradU.patchInternalField();
                    const tensorField pgradURight = pGradU.patchNeighbourField();

                    const vectorField pgradTLeft = pGradT.patchInternalField();
                    const vectorField pgradTRight = pGradT.patchNeighbourField();

                    // Geometry: call the raw cell-to-face vector by calling
                    // the base patch (cell-to-face) delta coefficient
                    vectorField pDeltaRLeft;
                    vectorField pDdeltaRRight;
                    if (U_.boundaryField()[patchi].type() == "cyclic")
                    {
                        // Work out the right delta from the cell-to-cell delta
                        // across the coupled patch and left delta
                        pDeltaRLeft = curPatch.fvPatch::delta();
                        pDdeltaRRight = pDeltaRLeft - curPatch.delta();
                    }
                    else
                    {
                        const vectorField faceCenter = pp.patch().Cf();
                        const vectorField pCellCenterLeft = pCellCenter.patchInternalField();
                        const vectorField pCellCenterRight = pCellCenter.patchNeighbourField();
                        pDeltaRLeft = faceCenter - pCellCenterLeft;
                        pDdeltaRRight = faceCenter - pCellCenterRight;
                    }

                    forAll(pp, facei)
                    {
                        Flux::evaluateFlux(
                            pRhoFlux[facei],
                            pRhoUFlux[facei],
                            pRhoEFlux[facei],

                            ppLeft[facei] + (pDeltaRLeft[facei] & pgradPLeft[facei]),

                            ppRight[facei] + (pDdeltaRRight[facei] & pgradPRight[facei]),

                            pULeft[facei] + (pDeltaRLeft[facei] & pgradULeft[facei]),

                            pURight[facei] + (pDdeltaRRight[facei] & pgradURight[facei]),

                            pTLeft[facei] + (pDeltaRLeft[facei] & pgradTLeft[facei]),

                            pTRight[facei] + (pDdeltaRRight[facei] & pgradTRight[facei]),

                            pR[facei],
                            pR[facei],
                            pCv[facei],
                            pCv[facei],
                            pSf[facei],
                            pMagSf[facei],
                            pBuiEps[facei]);
                    }
                }
                else if (
                    isType<wallFvPatch>(p_.mesh().boundary()[patchi]))
                {
                    forAll(pp, facei)
                    {
                        Flux::evaluateFlux(
                            pRhoFlux[facei],
                            pRhoUFlux[facei],
                            pRhoEFlux[facei],
                            pp[facei],
                            pp[facei],
                            pU[facei],
                            pU[facei],
                            pT[facei],
                            pT[facei],
                            pR[facei],
                            pR[facei],
                            pCv[facei],
                            pCv[facei],
                            pSf[facei],
                            pMagSf[facei],
                            1.0);
                    }
                }
                else
                {
                    const scalarField ppLeft =
                        p_.boundaryField()[patchi].patchInternalField();

                    const vectorField pULeft =
                        U_.boundaryField()[patchi].patchInternalField();

                    const scalarField pTLeft =
                        T_.boundaryField()[patchi].patchInternalField();

                    forAll(pp, facei)
                    {
                        Flux::evaluateFreestreamFlux(
                            pRhoFlux[facei],
                            pRhoUFlux[facei],
                            pRhoEFlux[facei],
                            ppLeft[facei],
                            pp[facei],
                            pULeft[facei],
                            pU[facei],
                            pTLeft[facei],
                            pT[facei],
                            pR[facei],
                            pR[facei],
                            pCv[facei],
                            pCv[facei],
                            pSf[facei],
                            pMagSf[facei],
                            1.0);
                    }
                }
            }
        }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
