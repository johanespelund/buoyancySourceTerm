/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Cg_ = coeffs().lookupOrDefault<scalar>("Cg", scalar(1)/0.85);
    Cphi_ = coeffs().lookupOrDefault<scalar>("Cphi", 0.3);
    THFM_ = coeffs().lookupOrDefault<word>("THFM", "SGDH");
    k_ = coeffs().lookupOrDefault<bool>("k", true);
    epsilon_ = coeffs().lookupOrDefault<bool>("epsilon", true);
    omega_ = coeffs().lookupOrDefault<bool>("omega", false);

    // Default v2_ to true when a v2-f model is detected
    v2_ = coeffs().lookupOrDefault<bool>("v2", isv2f_);

    Info<< "    Cg (1/Prt): " << Cg_ << " (" << 1/Cg_ << ")" << endl;
    Info<< "    THFM: " << THFM_ << endl;
    Info<< "    Applying source to the following fields:" << endl;
    Info<< "      k: " << k_ << endl;
    Info<< "      epsilon: " << epsilon_ << endl;
    Info<< "      omega: " << omega_ << endl;
    Info<< "      v2: " << v2_ << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::fv::buoyancyTurbSource::Gb(const volScalarField& rho) const
{
    const dimensionedScalar k0("k0", dimEnergy/dimMass, small);

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
            THFM_ == "SGDH"
          ?
            // SGDH formulation, Ref. [2]
            -turbulence_.nut()*Cg_*(g_ & fvc::grad(rho))
          :
            // GGDH formulation, Ref. [1]
            -1.5*Cg_*(turbulence_.nut()/(turbulence_.k() + k0))
           *((turbulence_.sigma() & fvc::grad(rho)) & g_)
        )
    );
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceEpsilon
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    const dimensionedScalar C1
    (
        "C1",
        dimless,
        turbulence_.coeffDict().lookupOrDefault<scalar>("C1", 1.44)
    );
    const volVectorField& U = turbulence_.U();
    const volScalarField& k = turbulence_.k();
    const dimensionedScalar k0(k.dimensions(), small);

    // Stratification angle factor C3 = tanh(|u_parallel / u_perp|), Ref. [3]
    const vector gHat(g_.value()/mag(g_.value()));

    const volScalarField uParallel(gHat & U);
    const volScalarField uPerp
    (
        mag(U - gHat*uParallel)
      + dimensionedScalar(dimVelocity, small)
    );

    const volScalarField GbField(Gb(rho));

    eqn -= fvm::SuSp(-C1*tanh(mag(uParallel/uPerp))*GbField/(k + k0), eqn.psi());
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceOmega
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    const volScalarField& nut = turbulence_.nut();
    const dimensionedScalar nut0(nut.dimensions(), small);
    const volScalarField& omega = eqn.psi();
    const dimensionedScalar omega0(omega.dimensions(), small);

    const scalar gamma
    (
        turbulence_.coeffDict().lookupOrDefault<scalar>("gamma", scalar(5)/9)
    );

    const volScalarField GbField(Gb(rho));

    eqn -= fvm::SuSp(-gamma/(nut + nut0)*GbField/(omega + omega0), omega);
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceK
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    const volScalarField& k = eqn.psi();
    const dimensionedScalar k0(k.dimensions(), small);
    const volScalarField GbField(Gb(rho));

    if (mesh().time().writeTime())
    {
        GbField.write();
    }

    eqn -= fvm::SuSp(-GbField/(k + k0), k);
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourcev2
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    const volScalarField& k = turbulence_.k();
    const dimensionedScalar k0(k.dimensions(), small);
    const volScalarField& v2 = eqn.psi();
    const volScalarField GbField(Gb(rho));

    // Buoyancy contribution to v2 scaled by the v2/k stress ratio,
    // using SuSp for implicit treatment under stable stratification (Gb < 0)
    eqn -= fvm::SuSp(-GbField/(k + k0), v2);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
    if (isv2f_ && v2_)
    {
        fields.append("v2");
    }

    return fields;
}


void Foam::fv::buoyancyTurbSource::addSup
(
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    const word& fieldName = field.name();

    if (fieldName == "k" && k_)
    {
        buoyancyTurbSourceK(rho, eqn);
    }
    else if (fieldName == "epsilon" && epsilon_)
    {
        buoyancyTurbSourceEpsilon(rho, eqn);
    }
    else if (fieldName == "omega" && omega_)
    {
        buoyancyTurbSourceOmega(rho, eqn);
    }
    else if (fieldName == "v2" && v2_)
    {
        buoyancyTurbSourcev2(rho, eqn);
    }
}


Foam::fv::buoyancyTurbSource::buoyancyTurbSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    g_(mesh.lookupObject<uniformDimensionedVectorField>("g")),
    turbulence_(mesh.lookupType<compressibleMomentumTransportModel>())
{
    // Detect if epsilon- or omega-based model is present
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

    // Detect v2-f model by checking for v2 and f fields
    isv2f_ =
        mesh.foundObject<volScalarField>("v2")
     && mesh.foundObject<volScalarField>("f");

    if (isv2f_)
    {
        Info<< "    Found v2 and f fields. Assuming v2-f model." << endl;
    }

    readCoeffs();
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

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
