/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "buoyancyTurbSource.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "fvcGrad.H"
#include "compressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(buoyancyTurbSource, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        buoyancyTurbSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::buoyancyTurbSource::readCoeffs()
{
    /*rhoName_ = coeffs().lookupOrDefault<word>("rho", "rho");*/
    /*alphatName_ = coeffs().lookupOrDefault<word>("alphat", "alphat");*/
    /*Tname_ = coeffs().lookupOrDefault<word>("T", "T");*/
    Cg_ = coeffs().lookupOrDefault<scalar>("Cg", 1/0.85);
    onlyApplyToK_ = coeffs().lookupOrDefault<bool>("onlyApplyToK", false);
    /*scalar beta_temp_ = coeffs().lookupOrDefault<scalar>("beta", -1);*/

    /*if (coeffs().found("beta"))*/
    /*{*/
    /*    beta_ = dimensionedScalar("beta", dimless/dimTemperature, beta_temp_);*/
    /*    Info << "  Using beta = " << beta_.value() << " for buoyancy source term." << nl;*/
    /*}*/
    /*else*/
    /*{*/
    /*    Info << "  No beta found. Using rho for buoyancy source term." << nl;*/
    /*    beta_ = dimensionedScalar("beta", dimless/dimTemperature, -1);*/
    /*}*/
}

Foam::tmp<Foam::volScalarField> Foam::fv::buoyancyTurbSource::B(const volScalarField&rho) const
{

    const auto& g = mesh().lookupObject<uniformDimensionedVectorField>("g");
    const auto& nut = turbulence_.nut();
    /*const volScalarField& k = turbulence_.k();*/
    /*auto k0 = dimensionedScalar("k0", k.dimensions(), SMALL);*/

    return tmp<Foam::volScalarField> 
    (
        new Foam::volScalarField
        (
            IOobject
            (
                "buoyancySource",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            /*nut() * Cg_ * (g & fvc::grad(rho) )/(k + k0)  // Computed buoyancy source expression*/
            nut() * Cg_ * (g & fvc::grad(rho) )   // Computed buoyancy source expression
        )
    );
}




// Apply source term to epsilon equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceEpsilon(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const dictionary& turbDict = turbulence_.coeffDict();
    const dimensionedScalar C1(turbDict.lookupOrDefault<scalar>("C1", 1.44));
    const volVectorField& U = turbulence_.U();
    const volScalarField& k = turbulence_.k();
    auto k0 = dimensionedScalar("k0", k.dimensions(), SMALL);

    // (BMA:Eq. 9)
    const vector gHat(g_.value()/mag(g_.value()));
 
    volScalarField v(gHat & U);
    volScalarField u
    (
        mag(U - gHat*v)
      + dimensionedScalar(dimVelocity, SMALL)
    );
 
    // (BMA:Eq. 6)
    const volScalarField _B = B(rho);

    eqn -= fvm::SuSp(C1*tanh(mag(v/u))*_B/(k + k0), eqn.psi());
}


// Apply source term to omega equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceOmega(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const volScalarField& nut = turbulence_.nut();
    /*const volScalarField& k = turbulence_.k();*/
    /*const volScalarField& omega = turbul*/
    const scalar gamma = 0.52;
    const volScalarField _B = B(rho);

    eqn -= gamma  / (nut + dimensionedScalar(nut.dimensions(), SMALL)) * _B;
    /*eqn -= fvm::SuSp(gamma  / (nut + dimensionedScalar(nut.dimensions(), SMALL)) * _B, eqn.psi());*/
}


// Apply source term to k equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceK(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const volScalarField& k = eqn.psi();
    const dimensionedScalar k0(k.dimensions(), SMALL);
    const volScalarField _B = B(rho);

    eqn -= fvm::SuSp(_B/(k + k0), k);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::buoyancyTurbSource::addSup(

                const volScalarField& rho,
  fvMatrix<scalar>& eqn, const word& fieldName) const
{
    if (fieldName == "k")
    {
        buoyancyTurbSourceK(rho, eqn);
    }
    else if (isEpsilon_ && fieldName == "epsilon")
    {
        buoyancyTurbSourceEpsilon(rho, eqn);
    }
    else if (!isEpsilon_ && fieldName == "omega")
    {
        buoyancyTurbSourceOmega(rho, eqn);
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported field name: " << fieldName << exit(FatalError);
    }
}

Foam::wordList Foam::fv::buoyancyTurbSource::addSupFields() const
{
    wordList fields;

    if (isEpsilon_ && !onlyApplyToK_)
    {
        fields.append("epsilon");
        Info << "   Adding buoyancy source to epsilon equation." << endl;
    }
    else if (!isEpsilon_ && !onlyApplyToK_)
    {
        fields.append("omega");
        Info << "   Adding buoyancy source to omega equation." << endl;
    }
    fields.append("k");

    return fields;
}



// Constructor with turbulence model loading and field detection
Foam::fv::buoyancyTurbSource::buoyancyTurbSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    /*beta_(dimensionedScalar("beta", dimless/dimTemperature, -1)),*/
    g_(mesh.lookupObject<uniformDimensionedVectorField>("g")),
    turbulence_(mesh.lookupType<momentumTransportModel>())
{
    readCoeffs();

    // Detect if epsilon or omega is present
    if (mesh.foundObject<volScalarField>("epsilon"))
    {
        isEpsilon_ = true;
    }
    else if (mesh.foundObject<volScalarField>("omega"))
    {
        isEpsilon_ = false;
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find an omega or epsilon field." << nl
            << "buoyancyTurbSource requires an omega- or epsilon-based model."
            << exit(FatalError);
    }
}


// Implementing pure virtual functions from fvModel as no-op or pass-through

bool Foam::fv::buoyancyTurbSource::movePoints()
{
    return true;
}


void Foam::fv::buoyancyTurbSource::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::buoyancyTurbSource::mapMesh(const polyMeshMap&)
{}


void Foam::fv::buoyancyTurbSource::distribute(const polyDistributionMap&)
{}


bool Foam::fv::buoyancyTurbSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //

