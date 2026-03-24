# buoyancySourceTerm

Adds a source term to transport equations in RANS models,
which accounts for preduction/destruction of turbulent kinetic energy
due to buoyancy effects.

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

