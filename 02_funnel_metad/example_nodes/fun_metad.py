#!/usr/bin/env python
# coding: utf-8

# # Funnel metadynamics
#
# This is a node to perform a funnel metadynamics simulation on a solvated protein-ligand complex.

# In[ ]:


import BioSimSpace as BSS


# Initialise the node and metadata.

# In[ ]:


node = BSS.Gateway.Node(
    "Perform funnel metadynamics on a solvated protein-ligand complex."
)
node.addAuthor(name="Dominykas Lukauskis", email="dominykas.lukauskis.19@ucl.ac.uk")


# Add the inputs to the node.

# In[ ]:


node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))

node.addInput(
    "runtime",
    BSS.Gateway.Time(
        help="The run time.",
        unit="nanoseconds",
        minimum=5 * BSS.Units.Time.nanosecond,
        maximum=2000 * BSS.Units.Time.nanosecond,
        default=1000 * BSS.Units.Time.nanosecond,
    ),
)

node.addInput(
    "hill_height",
    BSS.Gateway.Energy(
        help="The hill height.",
        unit="kj per mol",
        minimum=1 * BSS.Units.Energy.kj_per_mol,
        maximum=10 * BSS.Units.Energy.kj_per_mol,
        default=1.5 * BSS.Units.Energy.kj_per_mol,
    ),
)

node.addInput(
    "bias_factor",
    BSS.Gateway.Float(
        help="The bias factor for well-tempered metadynamics.",
        minimum=1.0,
        maximum=100.0,
        default=10.0,
    ),
)

node.addInput(
    "engine",
    BSS.Gateway.String(
        help="The molecular dynamics engine.",
        allowed=BSS.Metadynamics.engines(),
        default="OpenMM",
    ),
)

node.addInput(
    "work_dir",
    BSS.Gateway.String(
        help="The working directory for the simulation.", default="fun_metad"
    ),
)


# Add the outputs of the node.

# In[ ]:


node.addOutput(
    "final",
    BSS.Gateway.FileSet(
        help="The final system, in the original format plus a PDB file."
    ),
)


# Show the GUI so the user can set inputs interactively.

# In[ ]:


node.showControls()


# Create the system from the user specified input files.

# In[ ]:


system = BSS.IO.readMolecules(node.getInput("files"))


# Work out "best guess" funnel parameters based on the molecules in the system.

# In[ ]:


p0, p1 = BSS.Metadynamics.CollectiveVariable.makeFunnel(system)


# Create the funnel collective variable based on the funnel parameters.

# In[ ]:


cv = BSS.Metadynamics.CollectiveVariable.Funnel(
    p0, p1, upper_bound=BSS.Metadynamics.Bound(value=3.5 * BSS.Units.Length.nanometer)
)


# Create the metadynamics protocol, passing the collective variable and user-defined configuration parameters.

# In[ ]:


protocol = BSS.Protocol.Metadynamics(
    cv,
    runtime=node.getInput("runtime"),
    hill_height=node.getInput("hill_height"),
    bias_factor=node.getInput("bias_factor"),
    restart_interval=500000,
)


# Create a simulation process using the user-defined molecular-dynamics engine. This will automatically start the process.

# In[ ]:


engine = node.getInput("engine")
work_dir = node.getInput("work_dir")

if engine == "OpenMM":
    process = BSS.Process.OpenMM(system, protocol, platform="CUDA", work_dir=work_dir)
elif engine == "Amber":
    process = BSS.Process.Amber(
        system,
        protocol,
        exe="/home/model/MD-SOFTWARE/amber18-gnu-cu10/bin/pmemd.cuda",
        work_dir=work_dir,
    )
elif engine == "Gromacs":
    process = BSS.Process.Gromacs(system, protocol, work_dir=work_dir)

process.start()


# Wait for the process to finish and check for any errors.

# In[ ]:


process.wait()

# Print the last 10 lines of stdout and stderr if there was an error.
if process.isError():
    print(process.stdout(10))
    print(process.stderr(10))
    raise RuntimeError("The funnel metadynamics process exited with an error!")


# Get the system from the process, save to file, and bind to the output requirement.

# In[ ]:


final_system = process.getSystem()
formats = system.fileFormat() + ",pdb"
node.setOutput("final", BSS.IO.saveMolecules("final", final_system, formats))


# Validate the node.

# In[ ]:


node.validate()
