# Cell aggregation simulation

## DEMO
<img src="https://user-images.githubusercontent.com/40162543/45744016-63999480-bc38-11e8-9991-b57245957818.gif" width="50%" style="display:block;margin-left:auto;margin-right:auto;">

## Overview
Test code for the simulation of cell aggregation by chemotaxis.  Cells are represented as particles.  The motive force includes chemotaxis (proportinal to gradient of cAMP field) and Lennard Jones potential as volume exclusion.  The external field of cAMP is formed by point source secretion by each cells and diffusion.  Simulation is performed under periodic boundary condition.

## Requirements
- MATLAB (>= R2015b)
- Statistics and Machine Learning Toolbox

### How to run
Put cAMPsecretion_cellmigration_LJpotential.m file in MATLAB directory.  Execute `cAMPsecretion_cellmigration_LJpotential;`.  Initial cell position is random, and initial cAMP concentration is zero.  After the setting of initial condition, program will pause and require your input.  Press any key to start simulation.  If cells are too dense, error could occur due to too strong repulsive force by LJ potential.
