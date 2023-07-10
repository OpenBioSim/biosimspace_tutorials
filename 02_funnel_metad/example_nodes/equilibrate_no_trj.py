#!/usr/bin/env python
# coding: utf-8

# Author: Dominykas Lukauskis<br>
# Email:&nbsp;&nbsp; dominykas.lukauskis.19@ucl.ac.uk
#
# # Equilibration
#
# A node to perform equilibration of a molecular system. Saves the equilibration system to AMBER format files along with a DCD trajectory.

# In[ ]:


import BioSimSpace as BSS
import os

# set this just in case its missing in the env
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

# In[ ]:


node = BSS.Gateway.Node(
    "A node to perform equilibration. Saves the equlibrated molecular configuration to file."
)
node.addAuthor(
    name="Dominykas Lukauskis",
    email="dominykas.lukauskis.19@ucl.ac.uk",
    affiliation="University College London",
)
node.setLicense("GPLv3")


# Set the input requirements:

# In[ ]:


node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))

node.addInput(
    "engine",
    BSS.Gateway.String(
        help="The molecular dynamics engine.",
        allowed=BSS.MD.engines(),
        default="OpenMM",
    ),
)

node.addInput(
    "runtime",
    BSS.Gateway.Time(
        help="The run time.",
        unit="nanoseconds",
        minimum=0.02 * BSS.Units.Time.nanosecond,
        maximum=10 * BSS.Units.Time.nanosecond,
        default=2 * BSS.Units.Time.nanosecond,
    ),
)

node.addInput(
    "temperature_start",
    BSS.Gateway.Temperature(
        help="The initial temperature.",
        unit="kelvin",
        minimum=0 * BSS.Units.Temperature.kelvin,
        maximum=1000 * BSS.Units.Temperature.kelvin,
        default=0 * BSS.Units.Temperature.kelvin,
    ),
)

node.addInput(
    "temperature_end",
    BSS.Gateway.Temperature(
        help="The final temperature.",
        unit="kelvin",
        minimum=0 * BSS.Units.Temperature.kelvin,
        maximum=1000 * BSS.Units.Temperature.kelvin,
        default=300 * BSS.Units.Temperature.kelvin,
    ),
)

node.addInput(
    "restraint",
    BSS.Gateway.String(
        help="The type of restraint.",
        allowed=BSS.Protocol.Equilibration.restraints(),
        default="none",
    ),
)

# We now need to define the output of the node. In this case we will return a set of files representing the equilibrated molecular system in AMBER format and a single file containing the trajectory frames.

# In[ ]:


node.addOutput(
    "equilibrated", BSS.Gateway.FileSet(help="The equilibrated molecular system.")
)


# In[ ]:


node.showControls()


# Generate the molecular system.

# In[ ]:


system = BSS.IO.readMolecules(node.getInput("files"))


# Set up the equilibration protocol.
#
# (Note that the keyword arguments happen to have the same name as the input requirements. This need not be the case.)

# In[ ]:


protocol = BSS.Protocol.Equilibration(
    runtime=node.getInput("runtime"),
    temperature_start=node.getInput("temperature_start"),
    temperature_end=node.getInput("temperature_end"),
    restraint=node.getInput("restraint"),
)


# Start the molecular dynamics process.

# In[ ]:

engine = node.getInput("engine")

if engine == "OpenMM":
    process = BSS.Process.OpenMM(system, protocol, platform="CUDA")
elif engine == "Amber":
    process = BSS.Process.Amber(
        system, protocol, exe="/home/model/MD-SOFTWARE/amber18-gnu-cu10/bin/pmemd.cuda"
    )
elif engine == "Gromacs":
    process = BSS.Process.Gromacs(system, protocol)

process.start()


# We now wait for the process to finish, then check whether there were any errors before continuing. If errors were raised, then we raise an exception and print the last 10 lines of `stdout` and `stderr` to the user.

# In[ ]:


process.wait()

if process.isError():
    print(process.stdout(10))
    print(process.stdout(10))
    raise RuntimeError("The process exited with an error!")


# Get the equilibrated molecular system and write to file in the same format as the input, along with a PDB file to use as a topology file for the trajectory in [VMD](https://www.ks.uiuc.edu/Research/vmd/).

# In[ ]:


# Get the original file formats associated with the system.
# This is returned as a comma-separated string.
formats = system.fileFormat()

# Append PDB format, so we can visualise the trajectory with VMD.
formats += ",pdb"

# Write to file and bind to the output.
node.setOutput(
    "equilibrated", BSS.IO.saveMolecules("equilibrated", process.getSystem(), formats)
)


# Validate the node.

# In[ ]:


node.validate()
