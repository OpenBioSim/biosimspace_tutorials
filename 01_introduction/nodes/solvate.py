#!/usr/bin/env python
# coding: utf-8

# Author: Lester Hedges<br>
# Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk
#
# # Solvation
#
# A node to solvate a parameterised molecular system. An appropriate box size will be determine from the axis-aligned bounding box of the molecules in the system

# In[ ]:


import BioSimSpace as BSS


# In[ ]:


node = BSS.Gateway.Node(
    "A node to solvate a parameterised molecular system. An appropriate box size will be determine from the axis-aligned bounding box of the molecules in the system."
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
    "water_model",
    BSS.Gateway.String(
        help="The water model.", allowed=BSS.Solvent.waterModels(), default="tip3p"
    ),
)

node.addInput(
    "box",
    BSS.Gateway.String(
        help="The type of simulation box.", allowed=BSS.Box.boxTypes(), default="cubic"
    ),
)

node.addInput(
    "neutralise",
    BSS.Gateway.Boolean(
        help="Whether to neutralise the system. Ions will be added on top of any specified with 'ion_conc'",
        default=True,
    ),
)

node.addInput(
    "ion_conc",
    BSS.Gateway.Float(help="The ion concentration in mol per litre", default=0),
)


# We now need to define the output of the node. In this case we will return a set of files representing the parameterised molecule in the file format(s) of the original system.

# In[ ]:


node.addOutput("solvated", BSS.Gateway.FileSet(help="The solvated system."))


# In[ ]:


node.showControls()


# Load the molecular system.

# In[ ]:


system = BSS.IO.readMolecules(node.getInput("files"))


# Now work out the minimum box size for the molecules in the system.

# In[ ]:


# Get the minimium and maximum coordinates of the bounding box that
# minimally encloses the protein.
box_min, box_max = system.getAxisAlignedBoundingBox()

# Work out the box size from the difference in the coordinates.
box_size = [y - x for x, y in zip(box_min, box_max)]

# How much to pad each side of the molecule? (Nonbonded cutoff = 10 A)
padding = 15 * BSS.Units.Length.angstrom

# Work out an appropriate box. This will used in each dimension to ensure
# that the cutoff constraints are satisfied if the molecule rotates.
box_length = max(box_size) + 2 * padding


# Get the box parameters for the chosen box type.

# In[ ]:


box, angles = BSS.Box.generateBoxParameters(node.getInput("box"), box_length)


# Now solvate with the chosen water model and ion concentration.

# In[ ]:


solvated = BSS.Solvent.solvate(
    node.getInput("water_model"),
    molecule=system,
    box=box,
    angles=angles,
    is_neutral=node.getInput("neutralise"),
    ion_conc=node.getInput("ion_conc"),
)


# Write the solvated system to file in the same format as the original system.

# In[ ]:


node.setOutput(
    "solvated", BSS.IO.saveMolecules("solvated", solvated, system.fileFormat())
)


# Validate the node.

# In[ ]:


node.validate()
