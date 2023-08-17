/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "SIGMA.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
/*
template<class BasicTurbulenceModel>
tmp<volSymmTensorField> SIGMA<BasicTurbulenceModel>::Sd
(
    const volTensorField& gradU
) const
{
    return dev(symm(gradU & gradU));
}
*/

template<class BasicTurbulenceModel>
tmp<volScalarField> SIGMA<BasicTurbulenceModel>::k
(
    const volTensorField& gradU
) const
{
    volSymmTensorField G(symm(gradU.T()&gradU));

    volScalarField I1(tr(G));
    volScalarField I2(0.5*(sqr(I1)-tr(G&G)));
    volScalarField I3(det(G));

    volScalarField alpha1(sqr(I1)/9.0-I2/3.0);
    volScalarField alpha2(pow3(I1)/27.0-I1*I2/6.0+I3/2.0);
    volScalarField alpha3
    (
        acos(max(min(alpha2/pow(max(alpha1,sqr(small_)),1.5),1.0),-1.0))/3.0
    );

    scalar pi(constant::mathematical::pi);

    volScalarField sigma1
    (
        sqrt(max(I1/3.0+2.0*sqrt(alpha1)*cos(alpha3),small_))
    );
    volScalarField sigma2
    (
        sqrt(max(I1/3.0-2.0*sqrt(alpha1)*cos(pi/3.0+alpha3),zero_))
    );
    volScalarField sigma3
    (
        sqrt(max(I1/3.0-2.0*sqrt(alpha1)*cos(pi/3.0-alpha3),zero_))
    );

    volScalarField Dsigma(sigma3*(sigma1-sigma2)*(sigma2-sigma3)/sqr(sigma1));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            sqr(sqr(Csigma_)*this->delta()/Ck_)*sqr(Dsigma)
        )
    );
}


template<class BasicTurbulenceModel>
void SIGMA<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Ck_*this->delta()*sqrt(this->k(fvc::grad(this->U_)));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SIGMA<BasicTurbulenceModel>::SIGMA
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    ),

    Csigma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Csigma",
            this->coeffDict_,
            1.35
        )
    ),

    small_(dimensionedScalar("small", dimless/sqr(dimTime), SMALL)),
    zero_(dimensionedScalar("zero", dimless/sqr(dimTime), 0.0))
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SIGMA<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
        Csigma_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SIGMA<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}


template<class BasicTurbulenceModel>
void SIGMA<BasicTurbulenceModel>::correct()
{
    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
