#!/usr/bin/env python
# coding: utf-8

# Author: Lester Hedges<br>
# Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk
#
# # Parameterisation
#
# A node to perform parameterisation of a molecule loaded from PDB file. Saves the parameterised molecule in to AMBER format files.

# In[ ]:


import BioSimSpace as BSS


# In[ ]:


node = BSS.Gateway.Node(
    "A node to perform parameterisation of a molecule loaded from PDB file. Saves the parameterised molecule in to AMBER format files."
)
node.addAuthor(
    name="Lester Hedges",
    email="lester.hedges@bristol.ac.uk",
    affiliation="University of Bristol",
)
node.setLicense("GPLv3")


# Set the input requirements:

# In[ ]:


node.addInput(
    "pdb",
    BSS.Gateway.File(
        help="A Protein Data Bank (PDB) file containing a single molecule."
    ),
)

node.addInput(
    "forcefield",
    BSS.Gateway.String(
        help="The force field to parameterise the molecule with.",
        allowed=BSS.Parameters.forceFields(),
    ),
)

node.addInput(
    "water_model",
    BSS.Gateway.String(
        help="The water model to use for ion parameters.",
        allowed=BSS.Solvent.waterModels(),
        default="tip3p",
    ),
)

node.addInput(
    "pdb4amber",
    BSS.Gateway.Boolean(
        help="Whether to pre-process the PDB file using pdb4amber.", default=False
    ),
)


# We now need to define the output of the node. In this case we will return a set of files representing the parameterised molecule in AMBER format.

# In[ ]:


node.addOutput("parameterised", BSS.Gateway.FileSet(help="The parameterised molecule."))


# In[ ]:


node.showControls()


# Load the PDB file and pre-process with `pdb4amber` if requested. Since we assume a single molecule PDB, take the first molecule in the file.

# In[ ]:


molecule = BSS.IO.readPDB(node.getInput("pdb"), pdb4amber=node.getInput("pdb4amber"))[0]


# Perform the parameterisation using the chosen force field and ion water model. Note that we call the generic `BSS.Parameters.parameterise` function so that we can pass the force field name as an argument.

# In[ ]:


molecule = BSS.Parameters.parameterise(
    molecule, node.getInput("forcefield"), water_model=node.getInput("water_model")
).getMolecule()


# Now save the parameterise molecule to AMBER format files.

# In[ ]:


node.setOutput(
    "parameterised", BSS.IO.saveMolecules("parameterised", molecule, ["prm7", "rst7"])
)


# Validate the node.

# In[ ]:


node.validate()
