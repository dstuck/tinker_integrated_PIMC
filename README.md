# tinker_integrated_PIMC
PIMC code to calculate ZPE using force fields from tinker or other potentials

The code uses a harmonic oscillator propogator rather than the traditional free energy propogator to correct for anharmonic zero-point energy

A major caveat to this readme is that it is being put together 4 years after I have looked at or run the code

## Compiling
Install armadillo
- `brew install armadillo`
- set the environment variable `ARMADILLO_BASE` to the install location, i.e. '/usr/local/Cellar/armadillo/9.300.2'

Build with Make
- `make all`

## Running
Run the executable with a .pin input file, output location, and prm file for forcefields
- pimcTinker infile.pin outfile.out prmfile.prm

### Input format

There is a $molecule section which should mirror qchem inputs

The $params section is currently undocumented and needs to be inferred from Simulation.cpp. Importantly this section sets 
up variables about the system and defines which potential to use.

The $modes section allows the passing of normal modes as output by qchem (or at least qchem of 4 years ago)
