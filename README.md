# buoyancySourceTerm

Adds a source term to transport equations in RANS models,
which accounts for production/destruction of turbulent kinetic energy
due to buoyancy effects.

The user can choose between SGDH (simple gradient diffusion hypothesis)
and GGDH (generalised gradient diffusion hypothesis) for the turbulent
heat flux model. The term is always added to the `k` equation. It can
optionally be added to the `epsilon`, `omega`, and `v2` equations as well.

## Installation
Make sure OpenFOAM v12 is installed and sourced.
Create the installation directory, e.g.
```
mkdir -p $FOAM_RUN/../src
cd $FOAM_RUN/../src
```
Clone repository and compile
```
git clone https://github.com/johanespelund/buoyancySourceTerm.git
cd buoyancySourceTerm
wmake
```
## Usage
To use this fvModel, make sure that the following entry is in `system/controlDict`
```
libs
(
  "libbuoyancyTurbSourceFvModels.so"
);
```

Then add the model to `system/fvModels`. An example with all available
options and their default values is shown below:

```
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  12
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      fvModels;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

buoyancyTurbSource
{
    type        buoyancyTurbSource;

    // Turbulent heat flux model (THFM): SGDH or GGDH (default: SGDH)
    //   SGDH: simple gradient diffusion hypothesis
    //   GGDH: generalised gradient diffusion hypothesis
    THFM        SGDH;

    // Inverse turbulent Prandtl number (default: 1/0.85 ≈ 1.176)
    Cg          1.176;

    // Coefficient for buoyancy production term in GGDH (default: 0.3)
    Cphi        0.3;

    // Add source term to the k equation (default: true)
    k           true;

    // Add source term to the epsilon equation (default: true)
    epsilon     true;

    // Add source term to the omega equation (default: false)
    omega       false;

    // Apply stratification angle factor C3 = tanh(|u_parallel/u_perp|)
    // in the epsilon source term (default: true)
    tanhLimiter true;

    // Add source term to the v2 equation, v2-f models only
    // (default: auto-detected; true if v2 and f fields are present)
    // v2       true;
}

// ************************************************************************* //
```

