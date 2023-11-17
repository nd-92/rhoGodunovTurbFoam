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

\*---------------------------------------------------------------------------*/

#include "numericFlux.H"
#include "MDLimiter.H"
#include "cyclicFvPatch.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template <class Flux, class Limiter>
Foam::numericFlux<Flux, Limiter>::numericFlux(
    const volScalarField &p,
    const volVectorField &U,
    const volScalarField &T,
    // const surfaceScalarField &buiEps,
    const volScalarField &upwindingFactor,
    const basicThermo &thermo)
    : numericFluxBase<Flux>(p.mesh()),
      p_(p),
      U_(U),
      T_(T),
      buiEps_(fvc::interpolate(upwindingFactor)),
      thermo_(thermo),
      rhoFlux_(
          IOobject(
              "phi",
              this->mesh().time().timeName(),
              this->mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          (linearInterpolate(thermo_.rho() * U_) & this->mesh().Sf())),
      rhoUFlux_(
          IOobject(
              "rhoUFlux",
              this->mesh().time().timeName(),
              this->mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          rhoFlux_ * linearInterpolate(U_)),
      rhoEFlux_(
          IOobject(
              "rhoEFlux",
              this->mesh().time().timeName(),
              this->mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          rhoFlux_ * linearInterpolate(thermo.Cv() * T_ + 0.5 * magSqr(U_))),
      gradP(fvc::grad(p_)),
      gradU(fvc::grad(U_)),
      gradT(fvc::grad(T_))
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::computeFlux()
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

// ************************************************************************* //
