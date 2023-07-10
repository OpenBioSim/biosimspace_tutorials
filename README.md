# BioSimSpace tutorials

An introduction to BioSimSpace and advanced tutorials covering example use cases.
These tutorials have been tested with BioSimSpace 2023.3.0 and require the following dependencies

* Gromacs (tested with 2023.1)
* AmberTools (tested with 23.3) 
* PLUMED (tested with 2.9.0)
* cinnabar (tested with 0.3.0)
* alchemlyb (tested with 2.1.0)

We recommend using mamba to install the dependencies. 

```
mamba create -n  bsstutorials
mamba install -n bsstutorials -c openbiosim -c conda-forge biosimspace=2023.3.0 gromacs=2023.1 ambertools=23.3 plumed=2.9.0 cinnabar=0.3.0 alchemlyb=2.1.0
mamba activate bsstutorials
```

* 01 - [Introduction](01_introduction)
* 02 - [Funnel metadynamics](02_funnel_metad)
* 03 - [Steered molecular dynamics and ensemble building](03_steered_md)
* 04 - [Free-energy perturbation](04_fep)
