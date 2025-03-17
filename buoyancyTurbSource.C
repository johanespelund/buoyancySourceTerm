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
    onlyApplyToK_ = coeffs().lookupOrDefault<bool>("onlyApplyToK", false);
    THFM_ = coeffs().lookupOrDefault<word>("THFM", "SGDH");

    Info << "   Cg / Prt: " << Cg_ << " / " << 1/Cg_ << endl;
    Info << "   onlyApplyToK: " << onlyApplyToK_ << endl;
    Info << "   THFM: " << THFM_ << endl;
}

Foam::tmp<Foam::volScalarField> Foam::fv::buoyancyTurbSource::B(const volScalarField&rho) const
{

    const auto& g = mesh().lookupObject<uniformDimensionedVectorField>("g");
    const auto& nut = turbulence_.nut();
    const auto& R = turbulence_.sigma();
    const volScalarField& k = turbulence_.k();
    const volScalarField& epsilon = turbulence_.epsilon();
    /* scalar epsilonMin_ = turbulence_.epsilonMin(); */
    /*auto k0 = dimensionedScalar("k0", k.dimensions(), SMALL);*/
    /* volSymmTensorField R = ((2.0/3.0)*I)*k - (nut)*dev(twoSymm(fvc::grad(turbulence_.U()))); */

    /* Info << "   Calculating buoyancy source term." << endl; */

    /* volScalarField G = 0.9*(k/epsilon)*(R() & g & fvc::grad(rho)); // /(epsilon + SMALL) // GGDH */
    /* volScalarField G = (3/2)*Cg_*(nut()/k)*(R() & g & fvc::grad(rho)); // /(epsilon + SMALL) // GGDH */

    /* Info << "   Calculated buoyancy source term." << endl; */


    // Create vol fields that can be written for G as well
    
    /* Info << "   Writing buoyancy source term." << endl; */
    /* Foam::volScalarField Gvector */
    /* ( */
    /*     IOobject */
    /*     ( */
    /*         "buoyancySourceVector", */
    /*         mesh().time().timeName(), */
    /*         mesh(), */
    /*         IOobject::NO_READ, */
    /*         IOobject::NO_WRITE */
    /*     ), */
    /*     0.9*(k/epsilon)*(R() & g & fvc::grad(rho)) // /(epsilon + SMALL) // GGDH */
    /* ); */

    /* if (mesh().time().writeTime()) */
    /* { */
    /*     Gvector.write(); */
    /* } */

    /* Info << "   Writing R tensor." << endl; */
    /* // Same for R */
    /* Foam::volSymmTensorField Rvector */
    /* ( */
    /*     IOobject */
    /*     ( */
    /*         "Rvector", */
    /*         mesh().time().timeName(), */
    /*         mesh(), */
    /*         IOobject::NO_READ, */
    /*         IOobject::NO_WRITE */
    /*     ), */
    /*     /1* ((2.0/3.0)*I)*k - (nut)*dev(twoSymm(fvc::grad(turbulence_.U()))) *1/ */
    /*     R() */
    /* ); */
    
    /* if (mesh().time().writeTime()) */
    /* { */
    /*     Rvector.write(); */
    /* } */
    
    /* Info << "Returning buoyancy source term." << endl; */

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
            THFM_ == "SGDH" ? 
            nut() * Cg_ * (g & fvc::grad(rho))  // SGDH
            :   
            (3/2)*Cg_*(nut()/k)*(R() & g & fvc::grad(rho)) // /(epsilon + SMALL) // GGDH
            /* 0.9*(k/epsilon)*(R() & g & fvc::grad(rho)) // /(epsilon + SMALL) // GGDH */
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

    // Write the field to disk
    /* _B().write(); */
    // but only when time = writeTime!
    if (mesh().time().writeTime())
    {
        _B().write();
    }

    eqn -= fvm::SuSp(_B/(k + k0), k);
}

// Apply source term to f equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourcef(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const volScalarField& k = eqn.psi();
    const dimensionedScalar k0(k.dimensions(), SMALL);
    const volScalarField _B = B(rho);
    const scalar C2 = 0.3;

    eqn += C2*(_B/pow(k + k0,2), k);
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
    else if (isv2f_ && fieldName == "f")
    {
        buoyancyTurbSourcef(rho, eqn);
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
    if (isv2f_ && !onlyApplyToK_)
    {
        fields.append("f");
        Info << "   Adding buoyancy source to v2-f model." << endl;
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
    // Check for v2 and f fields to set isv2f_.
    if (mesh.foundObject<volScalarField>("v2") &&
        mesh.foundObject<volScalarField>("f"))
    {
        isv2f_ = true;
        Info << "   Found v2 and f fields. Assuming v2-f model." << endl;
    }
    else
    {
        isv2f_ = false;
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

