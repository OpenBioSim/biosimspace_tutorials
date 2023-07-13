# BioSimSpace tutorials

A suite of tutorials that provide an introduction to BioSimSpace, and examples of advanced functionality 
covering different scientific use cases. 
These tutorials have been tested with BioSimSpace 2023.3.0 and require the following dependencies

* Gromacs (tested with 2023.1)
* AmberTools (tested with 23.3) 
* PLUMED (tested with 2.9.0)
* cinnabar (tested with 0.3.0)
* alchemlyb (tested with 2.1.0)

# Installation instructions

## Option 1) Local installation using conda

This route can be used to install the tutorials on a computer.

We recommend using mamba to install the dependencies. 

```
conda install -c conda-forge mamba
mamba create -n  bsstutorials
mamba install -n bsstutorials -c openbiosim -c conda-forge biosimspace=2023.3.0 gromacs=2023.1 ambertools=23.3 plumed=2.9.0 cinnabar=0.3.0 alchemlyb=2.1.0
mamba activate bsstutorials
```

Next clone the repository

```
git clone https://github.com/OpenBioSim/biosimspace_tutorials
```

## Option 2) Use OpenBioSim's Jupyter Hub server

[OpenBioSim](https://www.openbiosim.org/) provides access to a Jupyter Hub server with a preinstalled python environment to run the tutorials. Use this route to run the tutorials from a compatible web-browser. Access to the server require a GitHub account.

Go to [try.openbiosim.org](https://try.openbiosim.org/)

And follow the instructions in the file TUTORIALS.txt.

# Tutorials suite

* 01 - [Introduction](01_introduction)
* 02 - [Funnel metadynamics](02_funnel_metad)
* 03 - [Steered molecular dynamics and ensemble building](03_steered_md)
* 04 - [Alchemical free energy calculations](04_fep)
