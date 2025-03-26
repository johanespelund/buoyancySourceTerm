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
    v2_ = coeffs().lookupOrDefault<bool>("v2", true);
    epsilon_ = coeffs().lookupOrDefault<bool>("epsilon", true);
    omega_ = coeffs().lookupOrDefault<bool>("omega", false);
    f_ = coeffs().lookupOrDefault<bool>("f", false);
    

    Info << "    Cg (1/Prt): " << Cg_ << " (" << 1/Cg_ << ")" << endl;
    Info << "    THFM: " << THFM_ << endl;
    Info << "    Applying source to following fields:" << endl;
    Info << "      k: " << k_ << endl;
    Info << "      v2: " << v2_ << endl;
    Info << "      epsilon: " << epsilon_ << endl;
    Info << "      omega: " << omega_ << endl;
    Info << "      f: " << f_ << endl;

}

Foam::tmp<Foam::volScalarField> Foam::fv::buoyancyTurbSource::B(const volScalarField&rho) const
{

    const uniformDimensionedVectorField& g = mesh().lookupObject<uniformDimensionedVectorField>("g");
    const volScalarField& nut = turbulence_.nut();
    /* const volSymmTensorField& R = turbulence_.sigma(); */
    /* const volScalarField& k = turbulence_.k(); */
    const volScalarField& epsilon = turbulence_.epsilon();

    // SMALL value with epsilon dims
    const dimensionedScalar eps0 = dimensionedScalar("SMALL", epsilon.dimensions(), SMALL);

    /* Info << "B(): Dimensions:" << endl; */
    /* Info << "   rho: " << rho.dimensions() << endl; */
    /* Info << "   g: " << g.dimensions() << endl; */
    /* Info << "   nut: " << nut.dimensions() << endl; */
    /* Info << "   R: " << R.dimensions() << endl; */
    /* Info << "   k: " << k.dimensions() << endl; */
    /* Info << "   epsilon: " << epsilon.dimensions() << endl; */

    /* const volScalarField SGDH = nut * Cg_ * (g & fvc::grad(rho)); */
    /* const volScalarField GGDH = (3/2)*Cg_*(k/epsilon)*(g & (R & fvc::grad(rho))); */
    /* const volVectorField grad_rho = fvc::grad(rho); */
    /* const volScalarField k_by_epsilon = k/epsilon; */
    /* const volVectorField R_times_grad_rho = turbulence_.sigma() & fvc::grad(rho); */
    

    /* Info << "   SGDH: " << SGDH.dimensions() << endl; */
    /* Info << "   GGDH: " << GGDH.dimensions() << endl; */
    /* Info << "   grad_rho: " << grad_rho.dimensions() << endl; */
    /* Info << "   k_by_epsilon: " << k_by_epsilon.dimensions() << endl; */
    /* Info << "   R_times_grad_rho: " << R_times_grad_rho.dimensions() << endl; */


    return tmp<Foam::volScalarField> 
    (
        new Foam::volScalarField
        (
            IOobject
            (
                "buoyancySource",
                mesh().time().name(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            THFM_ == "SGDH" ? 
            // Formulation from [3]
            nut * Cg_ * (g & fvc::grad(rho))
            :   
            // Formulation from [2]
            Cphi_*Cg_*(turbulence_.k()/(turbulence_.epsilon() + eps0))*(g & (turbulence_.sigma() & fvc::grad(rho))) // /(epsilon + SMALL) // GGDH
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

    // Check dimensions
    /* Info << " Dimensions of eqn: " << eqn.dimensions() << endl; */
    /* Info << " Dimensions of _B: " << _B().dimensions() << endl; */

    eqn -= fvm::SuSp(C1*tanh(mag(v/u))*_B/(k + k0), eqn.psi());
}

Foam::tmp<Foam::volScalarField> Foam::fv::buoyancyTurbSource::TsLocal() const
{
    /* const volScalarField& k = v2fModelPtr_->k(); */
    /* const volScalarField& epsilon = v2fModelPtr_->epsilon(); */
    /* const volScalarField& nu = v2fModelPtr_->nu(); */
    return
    Foam::tmp<Foam::volScalarField>
    (
        max
        (
            turbulence_.k()/turbulence_.epsilon(),
            6.0*sqrt(turbulence_.nu()/turbulence_.epsilon())
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::fv::buoyancyTurbSource::LsLocal() const
{
    const dimensionedScalar CL("CL", dimless, 0.23);   // ensure matches v2fCoeffs
    const dimensionedScalar Ceta("Ceta", dimless, 70); // ensure matches v2fCoeffs

    return 
    Foam::tmp<Foam::volScalarField>
    (
    CL * max
    (
        pow(turbulence_.k(), 1.5)/turbulence_.epsilon(),
        Ceta * pow025(pow3(turbulence_.nu())/turbulence_.epsilon())
      )
    );
}

// Apply source term to epsilon equation for v2f model
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceEpsilon_v2f(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const dictionary& turbDict = turbulence_.coeffDict();
    const dimensionedScalar C1(turbDict.lookupOrDefault<scalar>("C1", 1.44));
    /* const volVectorField& U = turbulence_.U(); */
    /* const volScalarField& k = turbulence_.k(); */

    /* const volScalarField& v2 = v2fModelPtr_->v2(); */
    const volScalarField Ceps1
    (
        1.4*(1.0 + 0.05*min(sqrt(turbulence_.k()/v2fModelPtr_->v2()), scalar(100.0)))
    );

    const volScalarField _B = B(rho);
    const volScalarField Seps = _B/TsLocal() * Ceps1;

    eqn -= fvm::SuSp(Seps/eqn.psi(), eqn.psi());
}

// Apply source term to omega equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceOmega(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const volScalarField& nut = turbulence_.nut();
    const scalar gamma = 0.52;
    const volScalarField _B = B(rho);

    /* Info << "   Adding buoyancy source to omega equation." << endl; */
    eqn -= gamma  / (nut + dimensionedScalar(nut.dimensions(), SMALL)) * _B;
}


// Apply source term to k equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceK(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const volScalarField& k = eqn.psi();
    const dimensionedScalar k0(k.dimensions(), SMALL);
    const volScalarField _B = B(rho);
    

    if (mesh().time().writeTime())
    {
        _B().write();
    }
    eqn -= fvm::SuSp(_B/(k + k0), k);
}

// Apply source term to f equation
void Foam::fv::buoyancyTurbSource::buoyancyTurbSourcef(const volScalarField&rho, fvMatrix<scalar>& eqn) const
{
    const volScalarField& k = turbulence_.k();
    const volScalarField _B = B(rho);
    const scalar C2 = 0.3;

    /* Info << "   Adding buoyancy source to f equation." << endl; */

    eqn += C2*_B/k;
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
    else if (isv2f_ && fieldName == "f" && f_)
    {
        buoyancyTurbSourcef(rho, eqn);
    }
    else if (isv2f_ && fieldName == "epsilon" && epsilon_) 
    {
        buoyancyTurbSourceEpsilon_v2f(rho, eqn);
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
    if (isv2f_ && f_)
    {
      fields.append("f");
    }
    if (isv2f_ && epsilon_)
    {
      fields.append("epsilon");
    }
    if (isv2f_ && v2_)
    {
      fields.append("v2");
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

