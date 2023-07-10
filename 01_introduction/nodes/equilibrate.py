#!/usr/bin/env python
# coding: utf-8

# Author: Lester Hedges<br>
# Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk
#
# # Equilibration
#
# A node to perform equilibration of a molecular system. Saves the equilibration system to AMBER format files along with a DCD trajectory.

# In[ ]:


import BioSimSpace as BSS


# In[ ]:


node = BSS.Gateway.Node(
    "A node to perform equilibration. Saves the equlibrated molecular configuration and trajectory to file."
)
node.addAuthor(
    name="Lester Hedges",
    email="lester.hedges@bristol.ac.uk",
    affiliation="University of Bristol",
)
node.setLicense("GPLv3")


# Set the input requirements:

# In[ ]:


node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))

node.addInput(
    "runtime",
    BSS.Gateway.Time(
        help="The run time.",
        unit="nanoseconds",
        minimum=0.02 * BSS.Units.Time.nanosecond,
        maximum=10 * BSS.Units.Time.nanosecond,
        default=0.02 * BSS.Units.Time.nanosecond,
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

node.addInput(
    "report_interval",
    BSS.Gateway.Integer(
        help="The number of integration steps between reporting output.",
        minimum=100,
        maximum=10000,
        default=100,
    ),
)

node.addInput(
    "restart_interval",
    BSS.Gateway.Integer(
        help="The number of integration steps between saving trajectory frames.",
        minimum=100,
        maximum=10000,
        default=500,
    ),
)

node.addInput(
    "engine",
    BSS.Gateway.String(
        help="The molecular dynamics engine.", allowed=BSS.MD.engines(), default="auto"
    ),
)


# We now need to define the output of the node. In this case we will return a set of files representing the equilibrated molecular system in AMBER format and a single file containing the trajectory frames.

# In[ ]:


node.addOutput(
    "equilibrated", BSS.Gateway.FileSet(help="The equilibrated molecular system.")
)
node.addOutput("trajectory", BSS.Gateway.File(help="The trajectory file."))


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


process = BSS.MD.run(system, protocol, engine=node.getInput("engine"))


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


# Get the trajectory, convert to MDTraj format and save to a DCD file.

# In[ ]:


# Get the BioSimSpace trajectory wrapper.
trajectory = process.getTrajectory()

# Extract the MDTraj object.
traj = trajectory.getTrajectory(format="mdtraj")

# Write to file.
traj.save("equilibrated.dcd")

# Bind to the output requirement.
node.setOutput("trajectory", "equilibrated.dcd")


# Validate the node.

# In[ ]:


node.validate()
