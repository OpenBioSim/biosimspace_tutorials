Author: Lester Hedges<br>
Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk

___Jupyter Recap___:
* Press Shift+Enter to execute a cell and move to the cell below.
* Press Ctrl+Enter to execute a cell and remain in that cell.
* Run a shell command on the underlying operating system by prefixing the command with an exclamation mark, !
* Remember that the flow is in the order that you execute cells, which is not necessarily linear in the notebook. Keep track of the numbers in brackets to the left of the cell!


# Molecular setup

## Introduction

In this section we will learn how to use BioSimSpace to set up a molecular system ready for simulation. Starting from a molecular topology in the form of a [Protein Data Bank](https://www.rcsb.org/) format file, we will learn how to parameterise molecules using different molecular [force fields](https://en.wikipedia.org/wiki/Force_field_(chemistry)), then solvate them using various [water models](https://en.wikipedia.org/wiki/Water_model).

### A note regarding molecular input

The starting point for many simulations is a molecular topology in the form of a [PDB](https://www.rcsb.org/) file. This file contains information regarding the structure of the molecule (its constituent residues and atoms), the layout of atoms in space (in the form of 3D atomic coordinates), and sometimes additional molecular information such as the formal charge of each atom. What this file does not contain is information describing how the atoms in the molecule _interact_, i.e. what are the functional forms and parameters for the terms in the molecular potential. This file is then used as the input to a _parameterisation engine_, which typically matches the atoms and residues against templates in order to _parameterise_ the molecule with a chosen force field. As such, the accuracy of the original topology is of critical importance: Atoms and residues _must_ have the correct names, and the topology _must_ be complete, i.e. no missing atoms.

Unfortunately, many tools do a poor job in preparing PDB files, e.g. having quirks with their naming conventions, excluding certain atoms, etc. Since it is impossible to account for all such inconsistencies, which often takes detailed knowledge of the particular system and tool in question, BioSimSpace takes the approach that the original files used to create a starting moleular system should be properly formatted from the outset. We don't want to make guesses as to what the user intended, or leave them confused if unexpected behaviour occurs later down the line.

If pre-processing of the PDB file is required, then we recommend using one of the following third-party tools:

* [pdb4amber](https://ambermd.org/tutorials/basic/tutorial9/index.php)
* [PDBFixer](https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html)

When present, we do provide rudimentary support for `pdb4amber` via the `BioSimSpace.IO.reaadPDB` function, where passing the `pdb4amber=True` argument will pre-process the file with `pdb4amber` prior to creating a molecular system. However, we choose only to support the _default_ options, since many are experimental and have can have undesirable knock-on effects, e.g. using the `--add-missing-atoms` option strips all chain identifiers from the molecule.

## Parameterisation

The [BioSimSpace.Parameters](https://biosimspace.org/api/index_Parameters.html) package provides support for parameterising molecules using three different engines:

* [AmberTools](https://ambermd.org/AmberTools.php) (Using the `tLEaP` and `antechamber` programs.)
* [gmx pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html) (Used as a fall-back for certain AMBER force fields when AmberTools isn't present.)
* [openff-toolkit](https://github.com/openforcefield/openff-toolkit) (The toolkit of the [Open Force Field Initiative](https://openforcefield.org/).)



Let's check if we need to download input files for this tutorial 


```python
from get_tutorial import download
download()
```

Let's load BioSimSpace and see what force fields are available:


```python
import BioSimSpace as BSS

BSS.Parameters.forceFields()
```


    


    INFO:rdkit:Enabling RDKit 2023.03.2 jupyter extensions





    ['ff03',
     'ff99',
     'ff99SB',
     'ff99SBildn',
     'ff14SB',
     'gaff',
     'gaff2',
     'openff_unconstrained-2.0.0-rc.1',
     'openff_unconstrained-1.2.0',
     'openff_unconstrained-2.1.0',
     'openff_unconstrained-1.3.1',
     'openff_unconstrained-1.1.0',
     'openff_unconstrained-2.0.0-rc.2',
     'openff_unconstrained-1.3.1-alpha.1',
     'openff_unconstrained-2.0.0',
     'openff_unconstrained-1.0.0-RC1',
     'openff_unconstrained-1.1.1',
     'openff_unconstrained-1.0.0-RC2',
     'openff_unconstrained-1.0.1',
     'openff_unconstrained-2.1.0-rc.1',
     'openff_unconstrained-1.3.0',
     'openff_unconstrained-1.0.0',
     'openff_unconstrained-1.2.1']



    /home/julien/miniforge3/envs/bsstutorials/lib/python3.10/site-packages/openff/interchange/smirnoff/_create.py:145: UserWarning: Automatically up-converting BondHandler from version 0.3 to 0.4. Consider manually upgrading this BondHandler (or <Bonds> section in an OFFXML file) to 0.4 or newer. For more details, see https://openforcefield.github.io/standards/standards/smirnoff/#bonds.
      _upconvert_bondhandler(force_field["Bonds"])
    /home/julien/miniforge3/envs/bsstutorials/lib/python3.10/site-packages/openff/interchange/smirnoff/_create.py:145: UserWarning: Automatically up-converting BondHandler from version 0.3 to 0.4. Consider manually upgrading this BondHandler (or <Bonds> section in an OFFXML file) to 0.4 or newer. For more details, see https://openforcefield.github.io/standards/standards/smirnoff/#bonds.
      _upconvert_bondhandler(force_field["Bonds"])


The supported force fields fall into two categories:

1) AMBER force fields:


```python
BSS.Parameters.amberForceFields()
```




    ['ff03', 'ff99', 'ff99SB', 'ff99SBildn', 'ff14SB', 'gaff', 'gaff2']



N.B. We currently don't support force fields from `AmberTools20` that use CMAP corrections.

2) Open Force Fields:


```python
BSS.Parameters.openForceFields()
```




    ['openff_unconstrained-2.0.0-rc.1',
     'openff_unconstrained-1.2.0',
     'openff_unconstrained-2.1.0',
     'openff_unconstrained-1.3.1',
     'openff_unconstrained-1.1.0',
     'openff_unconstrained-2.0.0-rc.2',
     'openff_unconstrained-1.3.1-alpha.1',
     'openff_unconstrained-2.0.0',
     'openff_unconstrained-1.0.0-RC1',
     'openff_unconstrained-1.1.1',
     'openff_unconstrained-1.0.0-RC2',
     'openff_unconstrained-1.0.1',
     'openff_unconstrained-2.1.0-rc.1',
     'openff_unconstrained-1.3.0',
     'openff_unconstrained-1.0.0',
     'openff_unconstrained-1.2.1']



N.B. We currently don't support the default _constrained_ versions of the force fields, since we require conversion via an intermediate [ParmEd](https://github.com/ParmEd/ParmEd) topology that needs explicit bond parameters. If required, constraints can be added at a later stage. This will hopefully be resolved future releases when direct translation from Open Force Field to BioSimSpace data structures should be possible.

N.B. The available Open Force Fields are determined dynamically at import time, so the list above might be different depending on what version of the `openff-toolkit` you have installed.

Let's load a small molecule and parameterise it with several supported force fields.


```python
# Load a methanol molecule from a PDB file. Since there is only a single
# molecule, we take the first item from the system that was created.
molecule = BSS.IO.readMolecules("inputs/methanol.pdb")[0]
```

As mentioned above, this is just a bare molecule that only contains information pertaining to the topology. To see this, we can query the _properties_ of the underlying [Sire](https://github.com/michellab/Sire) object.


```python
molecule._sire_object.propertyKeys()
```




    ['alt_loc',
     'space',
     'fileformat',
     'beta_factor',
     'element',
     'occupancy',
     'coordinates',
     'formal_charge',
     'is_ter',
     'connectivity',
     'is_het',
     'insert_code',
     'time']



We'll now parameterise the molecule with the [General AMBER force field](http://ambermd.org/antechamber/gaff.html), commonly known as GAFF. Behind the scenes this will set up and run the [antechamber](http://ambermd.org/tutorials/basic/tutorial4b/) and [tLEaP](https://ambermd.org/tutorials/pengfei/index.htm) programs from the [AmberTools](https://ambermd.org/AmberTools.php) suite. Depending on the input, `antechamber` might call out to `sqm` to perform a quantum chemistry calculation in order to calculate charges. Since this can be time consuming for a large molecule, all of the BioSimSpace parameterisation functions return a handle to a background process so that you can continue work interactively while you want for the the parameterisation completes.


```python
process = BSS.Parameters.gaff(molecule)
```

When you're ready to get the molecule, just call `.getMolecule()` on the process which will block until the parameterisation is complete, following which it will return a new molecule with force field parameters, or raise an exception if something went wrong.


```python
gaff_molecule = process.getMolecule()
```

N.B. If something went wrong, it can be useful look at the intermediate files within `process.workDir()` to see what errors were reported by the various programs that were run. A `README.txt` file in this directory will also tell you exactly what commands were run, and in what order.

Since this was just a small molecule and parameterisation was quick, we could have just returned the molecule from the process immediately using:


```python
gaff_molecule = BSS.Parameters.gaff(molecule).getMolecule()
```

N.B. When returning immediately any intermediate files will be lost unless the `work_dir` parameter was used to specify a working directory for the process.

Once again, we can query the underlying Sire object to see what properties are associated with the molecule:


```python
gaff_molecule._sire_object.propertyKeys()
```




    ['occupancy',
     'mass',
     'alt_loc',
     'beta_factor',
     'connectivity',
     'fileformat',
     'time',
     'intrascale',
     'is_het',
     'treechain',
     'atomtype',
     'gb_radii',
     'bond',
     'improper',
     'coordinates',
     'charge',
     'space',
     'formal_charge',
     'insert_code',
     'ambertype',
     'gb_radius_set',
     'forcefield',
     'LJ',
     'dihedral',
     'gb_screening',
     'element',
     'is_ter',
     'angle']



In addition to the properties loaded from the original PDB file we now have properties that relate to the force field parameters, such as `bond`, `angle`, and `dihedral`.

Note that when calling `.getMolecule()` BioSimSpace copies any additional properties from the parameterised system (created by loading the final output from the parameterisation process) back into the a copy of the original molecule, such that the _original_ topology is _preserved_. For example, while the parameterisation process may have renamed atoms/residues, or reordered atoms, the naming and ordering in the returned molecule will match the original that was passed in. As mentioned earlier, we don't deal with situations where the parameterisation engine _adds_ atoms that were missing from the original topology. In this case the parameterisation would fail, since the new topology is inconsistent with the original.

Let us now parameterise the same molecule using one of the Open Force Fields:


```python
openff_molecule = BSS.Parameters.openff_unconstrained_1_0_0(molecule).getMolecule()
```

We can now loop over atoms in the two parameterised molecules and compare properties. For example, we can see that the atomic charges are the same:


```python
for atom0, atom1 in zip(gaff_molecule.getAtoms(), openff_molecule.getAtoms()):
    print(atom0.name(), atom0.charge(), atom1.charge())
```

    C 0.1167 |e| 0.1167 |e|
    H1 0.0287 |e| 0.0287 |e|
    H2 0.0287 |e| 0.0287 |e|
    H3 0.0287 |e| 0.0287 |e|
    OH -0.5988 |e| -0.5988 |e|
    HO 0.3960 |e| 0.3960 |e|


To compare specific terms in the force field we can query the properties of the underlying Sire objects:


```python
# Get the bond potentials generated by GAFF.
gaff_molecule._sire_object.property("bond").potentials()
```




    [TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(4)} : 316.7 [r - 1.4233]^2 ),
     TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(3)} : 330.6 [r - 1.0969]^2 ),
     TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(2)} : 330.6 [r - 1.0969]^2 ),
     TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(1)} : 330.6 [r - 1.0969]^2 ),
     TwoAtomFunction( {CGIdx(0),Index(4)} <-> {CGIdx(0),Index(5)} : 371.4 [r - 0.973]^2 )]




```python
# Get the bond potentials generated by OpenFF.
openff_molecule._sire_object.property("bond").potentials()
```




    [TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(4)} : 334.571 [r - 1.41429]^2 ),
     TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(3)} : 379.047 [r - 1.09289]^2 ),
     TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(2)} : 379.047 [r - 1.09289]^2 ),
     TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(1)} : 379.047 [r - 1.09289]^2 ),
     TwoAtomFunction( {CGIdx(0),Index(4)} <-> {CGIdx(0),Index(5)} : 560.292 [r - 0.970769]^2 )]



As well as being able to parameterise a molecule loaded from file, BioSimSpace can also use a [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) string as the molecule parameter for any force field function. This can be useful if you want to highlight particular steroechemistry, which might not be possible in the intermediate file formats that are used behind the scenes duing the parameterisation, e.g. PDB. When using SMILES there is no constraint that the parameterised topology matches that of the original molecule, since we have not yet created a fully-fledged BioSimSpace molecule at the point at which we invoke the parameterisation. As such, it is perfectly acceptable for the parameterisation to add hydrogen atoms.

With this in mind, the two parameterisations shown above could also have been performed as follows:


```python
gaff_molecule = BSS.Parameters.gaff("CO").getMolecule()
openff_molecule = BSS.Parameters.openff_unconstrained_1_0_0("CO").getMolecule()
```

Here hydrogen atoms have been added during the parameterisation. While the atom layout and naming is different those of the previous examples, the parameters for the equivalent atoms are the same, e.g.:


```python
for atom0, atom1 in zip(gaff_molecule.getAtoms(), openff_molecule.getAtoms()):
    print(atom0.name(), atom0.charge(), atom1.charge())
```

    C1 0.1167 |e| 0.1167 |e|
    O2 -0.5988 |e| 0.0287 |e|
    H3 0.0287 |e| 0.0287 |e|
    H4 0.0287 |e| 0.0287 |e|
    H5 0.0287 |e| -0.5988 |e|
    H6 0.3960 |e| 0.3960 |e|


In addition to small molecules, BioSimSpace provides support for parameterising proteins using force fields from `AmberTools`, such as [ff14SB](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00255). (We don't currently support [ff19SB](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00591) due to the presence of CMAP terms, which we don't yet support in our parsers.)

As an example:


```python
protein = BSS.IO.readMolecules("inputs/2JJC.pdb")[0]
protein = BSS.Parameters.ff14SB(protein).getMolecule()
```

When molecules contain bound ions it is necessary to choose a water model for the ion parameters. This can be achieved by passing the `water_model` argument to any parameterisation function, where the named water model must match one of those described in the __Solvation__ section below, e.g. `water_model="tip3p"`. (An exception will be raised if bound ions are detected and no water model is chosen.) An optional `leap_commands` argument allows you to pass additional directives to the `tLEaP` program called by an AMBER protein force field function. These commands are added after the defaults allowing you to load custom force field parameters, etc.

## Solvation

The next stage in setting up a system ready for simulation is to solvate the molecule(s) in a box of water. The [BioSimSpace.Solvent](https://biosimspace.org/api/index_Solvent.html) package provides support for solvating with a variety of water models:


```python
BSS.Solvent.waterModels()
```




    ['spc', 'spce', 'tip3p', 'tip4p', 'tip5p']



N.B. At present we only support solvating in water.

Solvation is performed using the [gmx solvate](https://manual.gromacs.org/documentation/2018/onlinehelp/gmx-solvate.html) package. We currently don't support other solvation engines, since they require parameterisation as a pre-requisite, or include parameterisation as part of the solvation process itself, i.e. you can't decouple the two stages. We will hopefully overcome this shortcoming in future releases, since other engines would enable support for [more realistic salt concentrations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6078207/) and improved water templates, i.e. less vapour bubbles at box edges, which need to be properly equlibrated.

BioSimSpace provides support for both orthorhombic and triclinic simulation boxes, where appropriate box magnitudes and angles can be obtained using the [BioSimSpace.Box](https://biosimspace.org/api/index_Box.html) package. To see what pre-generated box types are available:


```python
BSS.Box.boxTypes()
```




    ['cubic',
     'rhombicDodecahedronHexagon',
     'rhombicDodecahedronSquare',
     'truncatedOctahedron']



For example, to get box parameters for a truncated octahedral box with an image distance of 10 nanometers.


```python
box, angles = BSS.Box.truncatedOctahedron(10 * BSS.Units.Length.nanometer)
print(box, angles)
```

    [100.0000 A, 100.0000 A, 100.0000 A] [70.5288 degrees, 109.4712 degrees, 70.5288 degrees]


When choosing a box for solvation it is important to ensure that it large enough to hold the molecule(s) of interest. In addition, when adding ions to neutralise the system or reach a desired ionic strength, then the system must be large enough that the cut-off used by [gmx genion](https://manual.gromacs.org/archive/5.0/programs/gmx-genion.html) when computing electrostatics is not more than twice the shortest box dimension. This can lead to a slight self-consistency issue, since it's sometimes not possible to know the number of ions that will need to be added a priori. Even if that were possible, then the random insertion of ions by `gmx genion` can still lead to issues if the choice of location means that they don't all manage to fit within the available volume.

As such, sometimes a little trial-and-error is needed to find an appropriate box size for the system in question. A good rule of thumb is to obtain the [axis-aligned bounding box](https://en.wikipedia.org/wiki/Bounding_volume) for the molecule(s) and add an appropriate buffer to the largest box dimension. For example:


```python
# Get the minimium and maximum coordinates of the bounding box that
# minimally encloses the protein.
box_min, box_max = protein.getAxisAlignedBoundingBox()

# Work out the box size from the difference in the coordinates.
box_size = [y - x for x, y in zip(box_min, box_max)]

# How much to pad each side of the protein? (Nonbonded cutoff = 10 A)
padding = 15 * BSS.Units.Length.angstrom

# Work out an appropriate box. This will used in each dimension to ensure
# that the cutoff constraints are satisfied if the molecule rotates.
box_length = max(box_size) + 2 * padding
```

Armed with this information, we can now solvate our protein in an appropriately sized cubic box. Here we will use the TIP3P water model:


```python
solvated = BSS.Solvent.tip3p(molecule=protein, box=3 * [box_length])
```

N.B. The `molecule` argument is optional. If ommited, then a pure water box will be generated.

By default, BioSimSpace will add counter-ions to neutralise the system. To see what ions were added we can use the built in `search` functionality:


```python
# Search for all free ions. As a simple search, we look for all molecules
# that only contain a single atom.
search = solvated.search("not mols with atomidx 2")

# Print all ions and their charge.
for ion in search.atoms():
    print(f"element = {ion.element()}, charge = {ion.charge()}")
```

    element = Sodium (Na, 11), charge = 1.0000 |e|
    element = Sodium (Na, 11), charge = 1.0000 |e|
    element = Sodium (Na, 11), charge = 1.0000 |e|
    element = Sodium (Na, 11), charge = 1.0000 |e|
    element = Sodium (Na, 11), charge = 1.0000 |e|
    element = Sodium (Na, 11), charge = 1.0000 |e|
    element = Sodium (Na, 11), charge = 1.0000 |e|


He're we see that `gmx genion` added 7 sodium ions. To confirm that the system was indeed neutralised, we can check its charge, as well as the charge of the original protein:


```python
print(f"solvated = {solvated.charge()}, protein = {protein.charge()}")
```

    solvated = -1.2183e-07 |e|, protein = -7.0000 |e|


We now have a solvated system ready for simulation. Let's visualise it with [BioSimSpace.Notebook.View](https://biosimspace.org/api/generated/BioSimSpace.Notebook.View.html#BioSimSpace.Notebook.View):


```python
view = BSS.Notebook.View(solvated)
view.system()
```


![Visualisation of the solvated system.](https://github.com/openbiosim/BioSimSpaceTutorials/blob/main/01_introduction/assets/02_solvated.png)   


Finally, let's save the system to file in AMBER format.


```python
BSS.IO.saveMolecules("solvated", solvated, ["prm7", "rst7"])
```




    ['/home/julien/projects/tmp/BioSimSpaceTutorials/01_introduction/solvated.prm7',
     '/home/julien/projects/tmp/BioSimSpaceTutorials/01_introduction/solvated.rst7']


