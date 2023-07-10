Author: Lester Hedges<br>
Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk

# BioSimSpace

## Introduction

Welcome to this workshop on [BioSimSpace](https://biosimspace.org), an _interoperable_ Python framework for biomolecular simulation. In this introductory session you will learn:

* What are the key concepts behind BioSimSpace.
* How to set up molecular systems ready for simulation.
* How to configure and run a range of molecular dynamics protocols using different simulation engines.
* How to write interoperable worflow components and run them in a variety of ways.

## What is BioSimSpace?

As a computational chemist you are likely overwhelmed by the amount of different software packages that are available to you. Having choice is a good thing, but too much can become a burden. I'm sure you have all come across at least one of the following:

* I know how to solve the problem with package X but I want to use package Y.
* How can I share my script with a collaborator who doesn't use the same software stack?
* How can I take advantage of the best tool for the job for different parts of my workflow?
* How can I compare methodology / results between simulation engines?

Solving these problems is the core goal of BioSimSpace. The wealth of fantastic software in our community is a real asset but _interoperability_ is currently a problem. Since there is no point reinventing the wheel, BioSimSpace is not an attempt to produce yet another molecular simulation package that reproduces all of the functionality from existing programs. This would result in just another tool for you to learn, along with yet another set of standards and formats. Instead, BioSimSpace is essentially just a set of _shims_, or bits of _glue_, that connect together existing software packages, allowing you to interact with them using a consistent Python interface.

## Why BioSimSpace?

By using BioSimSpace you will be less reliant on the use of brittle scripts to connect different software packages together. BioSimSpace builds on top of existing and open Python tools within the biomolecular community, e.g. [RDKit](https://www.rdkit.org/), [OpenMM](http://openmm.org/), [Open Force Field](https://github.com/openforcefield/openff-toolkit). As such, you are able to leverage the power of other packages, with which you may already be familiar, and to mix-and-match functionality where required.

With BioSimSpace you will be able to:

* Write generic workflow components _once_ in a package-agnostic language.
* Run the same script from the command-line, Jupyter, or within a workflow engine.
* Use the most suitable package that is availabe on your computer.
* Continue using your favourite package X but be able to share scripts with your collaborator who prefers package Y.
* Be able to take advantage of new software packages and hardware resources as and when they become available.

## What can BioSimSpace do?

BioSimSpace provides a suite of packages with a range of different functionality.

At present:

* File conversion: [BioSimSpace.IO](https://biosimspace.org/api/index_IO.html)
* Parameterisation: [BioSimSpace.Parameters](https://biosimspace.org/api/index_Parameters.html)
* Solvation: [BioSimSpace.Solvent](https://biosimspace.org/api/index_Solvent.html)
* Molecular dynamics: [BioSimSpace.Protocol](https://biosimspace.org/api/index_Protocol.html), [BioSimSpace.Process](https://biosimspace.org/api/index_Process.html), [BioSimSpace.MD](https://biosimspace.org/api/index_MD.html)
* Free-energy perturbation: [BioSimSpace.Align](https://biosimspace.org/api/index_Align.html), [BioSimSpace.FreeEnergy](https://biosimspace.org/api/index_FreeEnergy.html)
* Metadynamics: [BioSimSpace.Metadynamics](https://biosimspace.org/api/index_Metadynamics.html)
* Trajectory handling: [BioSimSpace.Trajectory](https://biosimspace.org/api/index_Trajectory.html)
* Interactive visualisation: [BioSimSpace.Notebook](https://biosimspace.org/api/index_Notebook.html)
* Workflow components: [BioSimSpace.Gateway](https://biosimspace.org/api/index_Gateway.html)

## Key concepts

Before getting started it's worth spending a little time covering a few of the key concepts of BioSimSpace.


### File conversion

While, broadly speaking, the different molecular dynamics engines offer a similar range of features, their interfaces are quite different. This makes it hard to take expertise in one package and immediately apply it to another. At the heart of this problem is the incompatibility between the molecular file formats used by the different packages. While they all contain the same information, i.e. how atoms are laid out in space and how they interact with each other, the structure of the files is very different. In order to provide interoperability betwen packages we will need to be able to read and write many different file formats, and be able to interconvert between them too.



First let's download input files for this tutorial


```python
from get_tutorial import download
download()
```

    Downloading inputs.tar.bz2 from openbiosim.org ...
    Extracting compressed tarball ...


Let's import the BioSimSpace Python package and see what we can do. For convenience, we'll rename the package to BSS to save us typing:


```python
import BioSimSpace as BSS
```


    


    INFO:rdkit:Enabling RDKit 2023.03.2 jupyter extensions


To see what file formats are supported by BioSimSpace, execute the cell below.


```python
BSS.IO.fileFormats()
```




    ['dcd',
     'gro87',
     'grotop',
     'mol2',
     'pdb',
     'prm7',
     'psf',
     'rst',
     'rst7',
     'sdf',
     'traj',
     'trr',
     'xtc']



Note that these refer to specific file _formats_, rather than file _extensions_. BioSimSpace doesn't care about file extensions, it's the _contents_ of the file that's important. 

If you aren't familiar with a particular format, you can get more information as follows, e.g.:


```python
BSS.IO.formatInfo("GroTop")
```




    'Gromacs Topology format files.'



The `BSS.IO.readMolecules` function is used to read molecular information from file. We've provided some example input files for you in the `inputs` directory. Let's take a look at some of these.

The `ala.crd` and `ala.top` files define a solvated alanine dipeptide system in AMBER format. Execute the cell below to see part of the topology file:


```python
!head -n 20 inputs/ala.top
```

    %VERSION  VERSION_STAMP = V0001.000  DATE = 06/30/15  11:44:23                  
    %FLAG TITLE                                                                     
    %FORMAT(20a4)                                                                   
    ACE                                                                             
    %FLAG POINTERS                                                                  
    %FORMAT(10I8)                                                                   
        1912       9    1902       9      25      11      43      24       0       0
        2619     633       9      11      24      13      21      20      10       1
           0       0       0       0       0       0       0       1      10       0
           0
    %FLAG ATOM_NAME                                                                 
    %FORMAT(20a4)                                                                   
    HH31CH3 HH32HH33C   O   N   H   CA  HA  CB  HB1 HB2 HB3 C   O   N   H   CH3 HH31
    HH32HH33O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  
    O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  
    H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   
    H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  
    O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  
    H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   
    H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  O   H1  H2  


Let's now read the molecules from file. The `BSS.IO.readMolecules` function automatically [globs](https://en.wikipedia.org/wiki/Glob_(programming)) the passed string, so wildcard matching can be used to determine the files. (Here the asterisk matches any characters, i.e. we are reading _all_ files in the `inputs` directory with `ala` as the file prefix.)


```python
system = BSS.IO.readMolecules("inputs/ala.*")
```

N.B. We could have explictly specified each file using a list of strings.

Note that we don't have to specify anything about the file format. BioSimSpace actually reads the files with _all_ of its inbuilt parsers in parallel. If a parser fails to read the file then it is immediately rejected and we move on to the next. As long as all of the files that are read contain a consistent topology, then BioSimSpace will be able to read them.

To see what file formats are associated with the system, run:


```python
system.fileFormat()
```




    'PRM7,RST7'



As expected, BioSimpace has detected that these were AMBER format topology and coordinate files.

The molecules are now loaded into a `System` object. BioSimSpace objects are thin wrappers around the equivalent objects from [Sire](https://github.com/michellab/Sire). For those that are familiar with Sire, you can get always access the underlying object directly using `system._sire_object`. This _private_ member is hidden from the user. Sire provides an extremely powerful and flexible set of tools for molecular manipulation and editing. While BioSimSpace directly exposes only a small subset of this functionality to the user, the full power of Sire is always easily available when needed.

We can query how many molecules there are as follows:


```python
system.nMolecules()
```




    631



To see how many water molecules there are:


```python
system.nWaterMolecules()
```




    630



To search for all nitrogen atoms in residues named `ALA`:


```python
system.search("element N and resname ALA")
```




    <BioSimSpace.SearchResult: nResults=1>



Create a new system using from the first 10 and last molecule in the system:


```python
new_system = (system[:10] + system[-1]).toSystem()
```

N.B. Information from the topology file pertaining to the molecular force field are stored internally as computer algebra expressions, allowing us to mathematically interconvert between different representations of the terms when writing to a different format. If desired, you could even edit the system to create your own unique force field parameters, although this beyond the scope of this tutorial.

Now that we have a molecular system, let's write it back to disk in a different format. Unsurprisingly, this is done using the `BSS.IO.saveMolecules` function. Execute the cell below to write the system to GROMACS format coordinate and topology files.


```python
BSS.IO.saveMolecules("ala", system, ["Gro87", "GroTop"])
```




    ['/home/julien/projects/tmp/BioSimSpaceTutorials/01_introduction/ala.gro',
     '/home/julien/projects/tmp/BioSimSpaceTutorials/01_introduction/ala.top']



Run the cell below to examine the start of the GROMACS topology file.


```python
!head -n 20 ala.top
```

    ; Gromacs Topology File written by Sire
    ; File written 07/07/23  23:52:48
    [ defaults ]
    ; nbfunc      comb-rule       gen-pairs      fudgeLJ     fudgeQQ
      1           2               yes            0.5         0.833333
    
    [ atomtypes ]
    ; name      at.num        mass      charge   ptype       sigma     epsilon
         C           6   12.010700    0.000000       A    0.339967    0.359824
        CT           6   12.010700    0.000000       A    0.339967    0.457730
        CX           6   12.010700    0.000000       A    0.339967    0.457730
         H           1    1.007940    0.000000       A    0.106908    0.065689
        H1           1    1.007940    0.000000       A    0.247135    0.065689
        HC           1    1.007940    0.000000       A    0.264953    0.065689
        HW           1    1.007940    0.000000       A    0.000000    0.000000
         N           7   14.006700    0.000000       A    0.325000    0.711280
         O           8   15.999400    0.000000       A    0.295992    0.878640
        OW           8   15.999400    0.000000       A    0.315061    0.636386
    
    [ moleculetype ]


### Topology preservation

One of the other pitfalls of working with different molecular simulation engines is that they often have quirks regarding naming conventions and atom ordering. This means that what you get back from a given program might not be the same as what you put in, making it tricky to cross-reference parts of the system that are of specific interest.

N.B. Differerent naming conventions are one of the hardest problems with interoperability.

BioSimSpace tries to preserve the intial molecular topology during any interaction with external tools so that it can be used as a consistent reference. For example, while we might need to rename the water topology to match the conventions of a particular molecular dynamics engine, we always copy the updated coordinates from a simulation back into the original system so that the naming that the user gets back is unchanged.

Another common topology feature that can be lost is chain identifiers. For example, a [Protein Data Bank](https://www.rcsb.org/) file might contain labels for chains, but these would be lost during parameterisation with the [AmberTools](https://ambermd.org/AmberTools.php) suite since it has no internal concept of chains.

As an example, consider the following PDB file:


```python
system_pdb = BSS.IO.readMolecules("inputs/1jr5.pdb")
```

This is read as a single molecule containing two chains:


```python
print(f"mols = {system_pdb.nMolecules()}, chains = {system_pdb.nChains()}")
```

    mols = 1, chains = 2


The molecule in the system has a topology but no force field information, as we can see by querying the _properties_ of the underlying Sire object:


```python
system_pdb[0]._sire_object.propertyKeys()
```




    ['beta_factor',
     'time',
     'coordinates',
     'insert_code',
     'is_ter',
     'fileformat',
     'is_het',
     'alt_loc',
     'formal_charge',
     'connectivity',
     'element',
     'space',
     'occupancy']



After parameterising the molecule with an AMBER protein force field, we end up with the following system:


```python
system_amber = BSS.IO.readMolecules(["inputs/1jr5.crd", "inputs/1jr5.top"])
```

In contrast, this is read as a two molecules with zero chains:


```python
print(f"mols = {system_amber.nMolecules()}, chains = {system_amber.nChains()}")
```

    mols = 2, chains = 0


Molecules in the parameterised system contain a different set of properties, including those pertaining to the force field.


```python
system_amber[0]._sire_object.propertyKeys()
```




    ['gb_radius_set',
     'parameters',
     'improper',
     'atomtype',
     'dihedral',
     'coordinates',
     'angle',
     'LJ',
     'intrascale',
     'space',
     'time',
     'treechain',
     'gb_screening',
     'connectivity',
     'mass',
     'fileformat',
     'bond',
     'gb_radii',
     'charge',
     'ambertype',
     'element',
     'forcefield']



If we want to preserve the topology of the original molecule, yet add the updated properties from the molecules in the new system, then we can make it _compatible_. BioSimSpace provides functionality to do this for you:


```python
# Extract the original molecule.
new_mol = system_pdb[0]

# Make it compatible with the new system, i.e. add new properties while
# preserving the topology.
new_mol.makeCompatibleWith(system_amber)
```

N.B. This is done implicitly whenever calling any BioSimSpace parameterisation function.

Let's now check the number of chains in the new molecule, as well as the properties that are associated with it.


```python
print(f"chains = {new_mol.nChains()}")
new_mol._sire_object.propertyKeys()
```

    chains = 2





    ['gb_radius_set',
     'improper',
     'atomtype',
     'dihedral',
     'occupancy',
     'coordinates',
     'angle',
     'LJ',
     'space',
     'intrascale',
     'time',
     'is_het',
     'treechain',
     'gb_screening',
     'connectivity',
     'mass',
     'insert_code',
     'fileformat',
     'gb_radii',
     'bond',
     'is_ter',
     'alt_loc',
     'charge',
     'ambertype',
     'element',
     'beta_factor',
     'forcefield',
     'formal_charge']


