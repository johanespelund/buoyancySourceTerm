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
    Cg_ = coeffs().lookupOrDefault<scalar>("Cg", 1/0.85);
    Cphi_ = coeffs().lookupOrDefault<scalar>("Cphi", 0.3); // Default value from [4]
    THFM_ = coeffs().lookupOrDefault<word>("THFM", "SGDH");
    k_ = coeffs().lookupOrDefault<bool>("k", true); 
    epsilon_ = coeffs().lookupOrDefault<bool>("epsilon", true);
    omega_ = coeffs().lookupOrDefault<bool>("omega", false);
    

    Info << "    Cg (1/Prt): " << Cg_ << " (" << 1/Cg_ << ")" << endl;
    Info << "    THFM: " << THFM_ << endl;
    Info << "    Applying source to following fields:" << endl;
    Info << "      k: " << k_ << endl;
    Info << "      epsilon: " << epsilon_ << endl;
    Info << "      omega: " << omega_ << endl;

}

Foam::tmp<Foam::volScalarField> Foam::fv::buoyancyTurbSource::Gb(const volScalarField&rho) const
{

    const uniformDimensionedVectorField& g = mesh().lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedScalar k0 = dimensionedScalar("k0", dimEnergy/dimMass, small);

    return tmp<Foam::volScalarField> 
    (
        new Foam::volScalarField
        (
            IOobject
            (
                "Gb",
                mesh().time().name(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            THFM_ == "SGDH" ? 
            // Formulation from [3]
            -turbulence_.nut() * Cg_ * (g & fvc::grad(rho))
            :   
            // Formulation from [2]
            -(3/2)*Cg_*(turbulence_.nut()/(turbulence_.k() + k0))*((turbulence_.sigma() & fvc::grad(rho)) & g) 
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
    auto k0 = dimensionedScalar("k0", k.dimensions(), small);

    // (BMA:Eq. 9)
    const vector gHat(g_.value()/mag(g_.value()));
 
    volScalarField v(gHat & U);
    volScalarField u
    (
        mag(U - gHat*v)
      + dimensionedScalar(dimVelocity, small)
    );
 
    const volScalarField _Gb = Gb(rho);

    eqn -= fvm::SuSp(-C1*tanh(mag(v/u))*_Gb/(k + k0), eqn.psi());
    /* eqn -= fvm::SuSp(-C1*_Gb/(k + k0), eqn.psi()); */
}



// Apply source term to omega equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceOmega(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const volScalarField& nut = turbulence_.nut();
    const scalar gamma = 0.52;
    const volScalarField _Gb = Gb(rho);

    eqn += gamma  / (nut + dimensionedScalar(nut.dimensions(), small)) * _Gb;
}


// Apply source term to k equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceK(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const volScalarField& k = eqn.psi();
    const dimensionedScalar k0(k.dimensions(), small);
    const volScalarField _Gb = Gb(rho);

    if (mesh().time().writeTime())
    {
        _Gb().write();
    }

    // Add source term using automatic treatment of explicit/implicit formulation
    // Ref. discussion here: 
    // https://www.cfd-online.com/Forums/openfoam-programming-development/182107-fvmatrix-fvoptions-susp-automatic-implicit-explicit-source-term-treatment.html

    eqn += -fvm::SuSp(-_Gb/(k + k0), k);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::buoyancyTurbSource::addSup(
                const volScalarField& rho,
  const volScalarField& field,
  fvMatrix<scalar>& eqn
) const
{

    word fieldName = field.name();
    if (fieldName == "k" && k_)
    {
        buoyancyTurbSourceK(rho, eqn);
    }
    else if (!isv2f_ && fieldName == "epsilon" && epsilon_)
    {
        buoyancyTurbSourceEpsilon(rho, eqn);
    }
    else if (fieldName == "omega" && omega_)
    {
        buoyancyTurbSourceOmega(rho, eqn);
    }
    else if (isv2f_) 
    {
        FatalErrorInFunction
            << "buoyancyTurbSource can not be used with v2f." << exit(FatalError);
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

    if (k_)
    {
      fields.append("k");
    }
    if (isEpsilon_ && epsilon_)
    {
      fields.append("epsilon");
    }
    if (!isEpsilon_ && omega_)
    {
      fields.append("omega");
    }

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
    turbulence_(mesh.lookupType<compressibleMomentumTransportModel>())
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
    // Check for v2 and f fields to set isv2f_.
    if (mesh.foundObject<volScalarField>("v2") &&
        mesh.foundObject<volScalarField>("f"))
    {
        isv2f_ = true;
        using Foam::RASModels::v2f;  // ensure namespace is available
        v2fModelPtr_ =
            dynamic_cast<const v2f<compressibleMomentumTransportModel>*>(&turbulence_);
        Info << "   Found v2 and f fields. Assuming v2-f model." << endl;

    }
    else
    {
        v2fModelPtr_ = nullptr;
        isv2f_ = false;
    }
    read(dict);
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

