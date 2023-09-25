# Funnel metadynamics

Author: Dominykas Lukauskis

Email: dominykas.lukauskis.19@ucl.ac.uk

Requirements:
* BioSimSpace
* OpenMM or Gromacs/Amber patched with PLUMED

This workshop tutorial walks you through how to setup funnel metadynamics (fun-metaD) simulations to estimate the absolute binding free energy (ABFE), how to visualise the funnel restraints and how to analyse the simulations using [BioSimSpace](https://biosimspace.org/index.html).

The workshop is divided into two notebooks. 

## 1. Simulation setup and restraint visualisation
For a step-by-step guide on how to visualise the funnel restraints and how to set up the fun-metaD simulations, follow [this Jupyter NB](01_bss-fun-metad-setup.ipynb).

## 2. Analysis
For a step-by-step guide on how to analyse fun-metaD simulations go [here](02_bss-fun-metad-analysis.ipynb). This notebook includes a discussion on how to think about convergence in fun-metaD simulations and how to get a robust ABFE estimate. 
