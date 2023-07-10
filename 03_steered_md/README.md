# Steered molecular dynamics

Author: Adele Hardie

Email: adele.hardie@ed.ac.uk

#### Requirements:
* BioSimSpace
* AMBER or GROMACS compiled with PLUMED
* An equilibrated starting system
* A target structure

This tutorial covers how to set up and run steered MD (sMD) simulations with BioSimSpace. In this example, the sMD trajectory is used to obtained seeds for equilibirum MD simulations for building a Markov State Model (MSM).

### 1. Set up and run sMD
An example of how to set up and run sMD with BioSimSpace is provided in [01_setup_sMD](01_setup_sMD.md). It covers some background for sMD, and creating steering protocols via BioSimSpace. Examples of both a single CV and a multi CV protocol are provided. A [notebook version](01_setup_sMD_AMBER.ipynb) exists as well.

### 2. Analysing trajectory and seeded MD
Once an sMD run is complete, the trajectory can be analysed and snapshots for seeded MD extracted [here](02_trajectory_analysis.ipynb). In this case we only look at the `COLVAR` file produced by PLUMED, but, depending on the system, other criteria might be important to decide whether the sMD run was successful. After the snapshots are extracted, they are used as starting points for equilibrium MD simulations, the basics of which are covered in the [introduction](../01_introduction).
