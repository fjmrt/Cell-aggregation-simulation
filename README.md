# Cell aggregation simulation

## DEMO
<img src="https://user-images.githubusercontent.com/40162543/45744016-63999480-bc38-11e8-9991-b57245957818.gif" width="50%" style="display:block;margin:auto;">

## Overview
Test code for the simulation of cell aggregation by chemotaxis.  Cells are represented as particles.  The motive force includes chemotactic force and Lennard-Jones potential as volume exclusion.  The direction of chemotactic force is gradient of cAMP field.  The strength of chemotactic force is constant, f0.  The dynamics is considered as over-damped condition, thus inertia term is ignored.  The external field of cAMP is formed by point source secretion by each cells and diffusion.  Simulation is performed under periodic boundary condition.

## Requirements
- MATLAB (>= R2015b)
- Statistics and Machine Learning Toolbox

## Parameters
- cellnum : number of cells.
- width : the length of edge of the simulated field.

### cAMP_secdegdif function
- k1 : cAMP production rate.
- k2 : cAMP degradation rate.
- DcAMP : diffusion constant of cAMP.

### cellmig function
- f0 : the strength of chemotactic force.
- epsilon, sigma : fitting parameters of [Lennerd-Jones potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential).

### How to run
Put cAMPsecretion_cellmigration_LJpotential.m file in MATLAB directory.  Execute `cAMPsecretion_cellmigration_LJpotential;`.  Initial cell position is random, and initial cAMP concentration is zero.  After the setting of initial condition, program will pause and require your input.  Press any key to start simulation.  If cells are too dense, error could occur due to too strong repulsive force by Lennard-Jones potential.
