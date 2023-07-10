# Part 2 - Setting up the system, visualising the funnel and preparing the fun-metaD simulations

In this part of the tutorial I will show you how to set up the BioSimSpace system, parameterising the protein and the ligand, as well as defining the simulation box, adding water and ions. Then, we will define a the funnel and visualise it using NGLview. Finally, we will setup the directories for the fun-metaD simulation itself.

## Part 2.1 - setup

First, let's download inputs for this tutorial


```python
from get_tutorial import download
download("01")
```

Now import BioSimSpace


```python
import BioSimSpace as BSS
import os
```

In order to parameterise the protein, you just need a tleap-ready protein PDB file. We will use Amber ff14sb forcefield.


```python
protein = BSS.IO.readMolecules("inputs_01/protein.pdb")
protein_water = protein.getWaterMolecules()

protein_search = protein.search("not water")
protein = protein_search.molecules().getResult(0)

protein = BSS.Parameters.ff14SB(protein, water_model="tip3p").getMolecule()
```

Parameterise the ligand using GAFF2.


```python
# parameterise the ligand
ligand = BSS.IO.readMolecules("inputs_01/ligand.mol2").getMolecule(0)
ligand = BSS.Parameters.gaff2(
    ligand, net_charge=BSS.Parameters.formalCharge(ligand)
).getMolecule()
```

Now we solvate the protein-ligand system using TIP3P water.


```python
# Find the dimensions of the protein
box_min, box_max = protein.getAxisAlignedBoundingBox()

# Work out the box size from the difference in the coordinates.
box_size = [y - x for x, y in zip(box_min, box_max)]

# how much to pad each side of the protein (nonbonded cutoff = 10 A)
padding = 15 * BSS.Units.Length.angstrom

box_length = max(box_size) + 2 * padding

cmplx = protein + ligand

solvated = BSS.Solvent.tip3p(molecule=cmplx.toSystem(), box=3 * [box_length])
```

Save the prepared structures.


```python
BSS.IO.saveMolecules("solvated", solvated, ["PDB", "RST7", "PRM7"])
```

## Part 2.2 - visualisation
Before we run any dynamics, let's visualise the funnel defined automatically by BSS.


```python
# Load the system.
system = BSS.IO.readMolecules(
    ["solvated.rst7", "solvated.prm7"]
)

# Create the funnel parameters.
p0, p1 = BSS.Metadynamics.CollectiveVariable.makeFunnel(system)

# Define the collective variable.
funnel_cv = BSS.Metadynamics.CollectiveVariable.Funnel(p0, p1)

# Create a view.
view = BSS.Metadynamics.CollectiveVariable.viewFunnel(system, funnel_cv)

view.molecules([1, 0, -1])
```

Lets check what the funnel looks with a larger radius.


```python
# Load the system.
system = BSS.IO.readMolecules(
    ["solvated.rst7", "solvated.prm7"]
)

# Create the funnel parameters.
p0, p1 = BSS.Metadynamics.CollectiveVariable.makeFunnel(system)

# Define the collective variable.
funnel_cv = BSS.Metadynamics.CollectiveVariable.Funnel(
    p0, p1, width=1.5 * BSS.Units.Length.nanometer
)

# Create a view.
view = BSS.Metadynamics.CollectiveVariable.viewFunnel(system, funnel_cv)

view.molecules([1, 0, -1])
```

In case the funnel was assigned incorrectly, you can open the structure file and use lists of atom IDs to define the p0 and p1 points.

Here I will load the same structure file with the p0 and p1 points defined by me.


```python
# Load the system.
system = BSS.IO.readMolecules(["solvated.pdb"])

# Create the funnel parameters.
p0, p1 = [2624], [584, 1307, 1474, 1887]

# Define the collective variable.
funnel_cv = BSS.Metadynamics.CollectiveVariable.Funnel(p0, p1)

# Create a view.
view = BSS.Metadynamics.CollectiveVariable.viewFunnel(system, funnel_cv)

view.molecules([1, 0, -1])
```

# Part 2.3 - Equilibrate the system with OpenMM
Load the prepared system first.


```python
solvated = BSS.IO.readMolecules(
    ["solvated.rst7", "solvated.prm7"]
)
```

The process below will take ~ 5 minutes if run on openbiosim's default Jupyter server. Running the notebook on a GPU powered computer should accelerate the calculations significantly


```python
protocol = BSS.Protocol.Minimisation()
process = BSS.Process.OpenMM(solvated, protocol)
# Start the process in the background.
process.start()

# Wait for the process to finish.
process.wait()
```

Get the minimised system and save the structure.


```python
minimised = process.getSystem()
# BSS.IO.saveMolecules('minimised',minimised, ['PDB','RST7','PRM7']) # if you want to save the output
```

Equilibrate (note the very short run time). This will take ca. 3 min on openbiosim's Jupyter server.


```python
protocol = BSS.Protocol.Equilibration(runtime=0.001 * BSS.Units.Time.nanosecond)
process = BSS.Process.OpenMM(minimised, protocol)

# Start the process in the background.
process.start()

# Wait for the process to finish.
process.wait()
```

Get the equilibrated system and save the structure


```python
equilibrated = process.getSystem()

# BSS.IO.saveMolecules('equilibrated',equilibrated, ['PDB','RST7','PRM7']) # if you want to save the output
```

## Part 2.4 - Setup and run fun-metaD

As I showed in the visualisation part, select the p0, p1 points for the funnel definition.


```python
p0, p1 = BSS.Metadynamics.CollectiveVariable.makeFunnel(equilibrated)
```

Write the funnel CV for the simulation. Additionally, we will set a lower upper_bound value.


```python
new_upper_bound = BSS.Metadynamics.Bound(value=3.5 * BSS.Units.Length.nanometer)

funnel_cv = BSS.Metadynamics.CollectiveVariable.Funnel(
    p0, p1, upper_bound=new_upper_bound
)
```

Write the metadynamics protocol.


```python
protocol = BSS.Protocol.Metadynamics(
    funnel_cv,
    runtime=0.01 * BSS.Units.Time.nanosecond,
    hill_height=1.5 * BSS.Units.Energy.kj_per_mol,
    restart_interval=1000,
    bias_factor=10,
)
```

Run the metadynamics simulation.


```python
process = BSS.Process.OpenMM(
    equilibrated, protocol, work_dir="fun-metaD-work_dir"
)

process.start()

process.wait()
```

Save the final structures.


```python
finished = process.getSystem()

# BSS.IO.saveMolecules('final', finished, ['PDB','RST7','PRM7'])  # if you want to save the output
```

## Closing thoughts

The purpose of this notebook has been to demonstrate the logic behind using BioSimSpace to setup and run your funnel metadynamics simulations, as well as visualising the suggested funnel restraints. 

In a real use case, you wouldn't do all of these steps in an interactive notebook. The funnel metadynamics simulations usually require 500-2000 ns of simulation time, which takes between 3 to 14 days for moderately sized systems. You want to submit a script that does the setting up, equilibration and production simulations all in one go. With Lester's help, I've prepared a set of nodes and a LSF submit script for doing all of the steps described above. Check out the `example_node` directory.
