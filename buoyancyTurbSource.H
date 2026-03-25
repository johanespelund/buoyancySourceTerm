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

Class
    Foam::fv::buoyancyTurbSource

Description
    Applies a buoyancy production term to the turbulence transport equations.
    The user can choose between SGDH (simple gradient diffusion hypothesis)
    and GGDH (generalised gradient diffusion hypothesis) for the turbulent
    heat flux model (THFM).

    The term is always added to the k equation. The user can choose to add
    the term to the epsilon, omega and v2 equations as well.

    For the epsilon equation, a stratification angle factor
    C3 = tanh(|u_parallel / u_perp|) (Ref. [3]) is applied by default.
    This can be disabled by setting \c tanhLimiter to \c false in the
    model dictionary.

    [1] Ince, N. Z., and B. E. Launder. 1989. "On the Computation of
        Buoyancy-Driven Turbulent Flows in Rectangular Enclosures."
        International Journal of Heat and Fluid Flow 10 (2): 110-17.
        https://doi.org/10.1016/0142-727X(89)90003-9.

    [2] Tieszen, S., A. Ooi, P. Durbin, and M. Behnia. 1998. "Modeling of
        Natural Convection Heat Transfer." Proceedings of the Summer Program.

    [3] Katsamis, C., T. Craft, H. Iacovides, and J. C. Uribe. 2022.
        "Use of 2-D and 3-D Unsteady RANS in the Computation of Wall
        Bounded Buoyant Flows." International Journal of Heat and Fluid
        Flow 93: 108914.
        https://doi.org/10.1016/j.ijheatfluidflow.2021.108914.

SourceFiles
    buoyancyTurbSource.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancyTurbSource_H
#define buoyancyTurbSource_H

#include "fvModel.H"
#include "compressibleMomentumTransportModel.H"
#include "uniformDimensionedFields.H"
#include "RASModel.H"

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class buoyancyTurbSource Declaration
\*---------------------------------------------------------------------------*/

class buoyancyTurbSource
:
    public fvModel
{
    // Private Data

        //- Flag indicating if an epsilon-based model is used
        bool isEpsilon_;

        //- Flag indicating if the v2-f model is used
        bool isv2f_;

        //- Flag: add source term to k equation
        bool k_;

        //- Flag: add source term to epsilon equation
        bool epsilon_;

        //- Flag: add source term to omega equation
        bool omega_;

        //- Flag: add source term to v2 equation (v2-f model only)
        bool v2_;

        //- Inverse turbulent Prandtl number
        scalar Cg_;

        //- Coefficient for the buoyancy production term in GGDH
        scalar Cphi_;

        //- Turbulent heat flux model (THFM) type.
        //  Either SGDH (simple gradient diffusion hypothesis, default) or
        //  GGDH (generalised gradient diffusion hypothesis).
        word THFM_;

        //- Flag: apply tanh limiter (stratification angle factor C3) in
        //  the epsilon source term
        bool tanhLimiter_;

        //- Gravitational acceleration [m/s^2]
        const uniformDimensionedVectorField& g_;

        //- Reference to the turbulence model
        const compressibleMomentumTransportModel& turbulence_;


    // Private Member Functions

        //- Read coefficients from dictionary
        void readCoeffs();

        //- Compute buoyancy production term [kg/(m s^3)]
        tmp<volScalarField> Gb(const volScalarField& rho) const;

        //- Apply buoyancy source to epsilon equation
        void buoyancyTurbSourceEpsilon
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn
        ) const;

        //- Apply buoyancy source to omega equation
        void buoyancyTurbSourceOmega
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn
        ) const;

        //- Apply buoyancy source to k equation
        void buoyancyTurbSourceK
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn
        ) const;

        //- Apply buoyancy source to v2 equation (v2-f model)
        void buoyancyTurbSourcev2
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn
        ) const;


public:

    //- Runtime type information
    TypeName("buoyancyTurbSource");

    // Constructors

        //- Construct from explicit source name and mesh
        buoyancyTurbSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        buoyancyTurbSource(const buoyancyTurbSource&) = delete;

    // Member Functions

        //- Return the list of fields for which the fvModel adds source term
        virtual wordList addSupFields() const;

        //- Add source term to turbulence transport equation
        virtual void addSup
        (
            const volScalarField& rho,
            const volScalarField& field,
            fvMatrix<scalar>& eqn
        ) const override;

        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints() override;

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&) override;

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&) override;

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&) override;

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict) override;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyancyTurbSource&) = delete;
};

} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
