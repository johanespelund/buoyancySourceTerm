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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::buoyancyTurbSource

Description
    Applies a source term to the turbulence model equations (k, epsilon, omega)
    to account for buoyancy effects in incompressible, compressible, and 
    multiphase flows.

SourceFiles
    buoyancyTurbSource.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancyTurbSource_H
#define buoyancyTurbSource_H

#include "fvModel.H"
#include "compressibleMomentumTransportModel.H"
#include "uniformDimensionedFields.H"

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

        //- Flag indicating if epsilon or omega model is used
        bool isEpsilon_;

        //- Flag indicatin if source term should only be added to k eqn.
        bool onlyApplyToK_; 


        //- Thermal expansion coefficient [1/K]
        /*dimensionedScalar beta_;*/

        //- Intervser turbulent Prandtl number
        scalar Cg_; 

        //- Gravitational acceleration vector [m/s^2]
        const uniformDimensionedVectorField& g_;

        //- Reference to the turbulence model
        const momentumTransportModel& turbulence_;

    // Private Member Functions

        //- Reads coefficients from dictionary
        void readCoeffs();

        //- Compute buoyancy production term for incompressible cases [m2/s3]
        tmp<volScalarField> B(const volScalarField& rho) const;
        /* volScalarField B() const; */

        //- Apply buoyancyTurbSource to epsilon field
        void buoyancyTurbSourceEpsilon(const volScalarField& rho, fvMatrix<scalar>& eqn) const;

        //- Apply buoyancyTurbSource to omega field
        void buoyancyTurbSourceOmega(const volScalarField& rho, fvMatrix<scalar>& eqn) const;

        //- Apply buoyancyTurbSource to k field
        void buoyancyTurbSourceK(const volScalarField& rho, fvMatrix<scalar>& eqn) const;

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

        //- Add explicit source term to k, epsilon, or omega equation
        /* virtual void addSup(fvMatrix<scalar>& eqn, const word& fieldName) const override; */
        virtual void addSup(
                const volScalarField& rho,
                fvMatrix<scalar>& eqn, 
                const word& fieldName) const override;

        //- Return the list of fields for which the fvModel adds source term
        virtual wordList addSupFields() const override;

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

