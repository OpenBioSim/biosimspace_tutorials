#!/usr/bin/env python
# coding: utf-8

# Author: Dominykas Lukaukis<br>
# Email:&nbsp;&nbsp; dominykas.lukauskis.19@ucl.ac.uk
#
# # Parameterisation
#
# A node to perform parameterisation of a molecule loaded from PDB file. Saves the parameterised molecule in to AMBER format files.

# In[ ]:


import BioSimSpace as BSS


# In[ ]:


node = BSS.Gateway.Node(
    "A node to perform parameterisation and solvation of protein and ligand molecules loaded from PDB/MOL2 files. Saves the parameterised and solvated system to AMBER format files."
)
node.addAuthor(
    name="Dominykas Lukauskis",
    email="dominykas.lukauskis.19@ucl.ac.uk",
    affiliation="University College London",
)
node.setLicense("GPLv3")


# Set the input requirements:

# In[ ]:


node.addInput(
    "protein_pdb",
    BSS.Gateway.File(
        help="A Protein Data Bank (PDB) file containing a single molecule."
    ),
)

node.addInput(
    "ligand_file",
    BSS.Gateway.File(
        help="A Protein Data Bank (PDB) or Tripos MOL2 file containing a single organic molecule."
    ),
)

node.addInput(
    "protein_forcefield",
    BSS.Gateway.String(
        help="The protein force field to parameterise the molecule with.",
        allowed=BSS.Parameters.forceFields(),
        default="ff14SB",
    ),
)

node.addInput(
    "ligand_forcefield",
    BSS.Gateway.String(
        help="The ligand force field to parameterise the molecule with.",
        allowed=BSS.Parameters.forceFields(),
        default="gaff2",
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
        help="Whether to pre-process the protein PDB file using pdb4amber.",
        default=False,
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


# We now need to define the output of the node. In this case we will return a set of files representing the parameterised molecule in AMBER format.

# In[ ]:


node.addOutput(
    "solvated", BSS.Gateway.FileSet(help="The parameterised and solvated system.")
)


# In[ ]:


node.showControls()


# Load the PDB file and pre-process with `pdb4amber` if requested. Since we assume a single molecule PDB, take the first molecule in the file.

# In[ ]:


protein = BSS.IO.readPDB(
    node.getInput("protein_pdb"), pdb4amber=node.getInput("pdb4amber")
)[0]


# Perform the parameterisation using the chosen force field and ion water model. Note that we call the generic `BSS.Parameters.parameterise` function so that we can pass the force field name as an argument.

# In[ ]:


protein = BSS.Parameters.parameterise(
    protein,
    node.getInput("protein_forcefield"),
    water_model=node.getInput("water_model"),
).getMolecule()

# Load the PDB file and pre-process with `pdb4amber` if requested. Since we assume a single molecule PDB, take the first molecule in the file.

# In[ ]:


ligand = BSS.IO.readPDB(node.getInput("ligand_file"))[0]


# Perform the parameterisation using the chosen force field and ion water model. Note that we call the generic `BSS.Parameters.parameterise` function so that we can pass the force field name as an argument.

# In[ ]:

############ net charge assignment? ##############
ligand = BSS.Parameters.parameterise(
    ligand,
    node.getInput("ligand_forcefield"),
).getMolecule()

# Combine the protein and the ligand

cmplx = protein + ligand

# Now work out the minimum box size for the molecules in the system.

# In[ ]:


# Get the minimium and maximum coordinates of the bounding box that
# minimally encloses the protein.
box_min, box_max = cmplx.getAxisAlignedBoundingBox()

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
    molecule=cmplx,
    box=box,
    angles=angles,
    is_neutral=node.getInput("neutralise"),
    ion_conc=node.getInput("ion_conc"),
)


# Write the solvated system to file in PDB and Amber format.

# In[ ]:


node.setOutput(
    "solvated", BSS.IO.saveMolecules("solvated", solvated, ["PDB", "RST7", "PRM7"])
)


# Validate the node.

# In[ ]:


node.validate()
