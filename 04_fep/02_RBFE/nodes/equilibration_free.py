#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import sys
from pathlib import Path

# In[ ]:


# Helper function
def runProcess(system, protocol, engine="GROMACS", AMBER_path=None, pmemd=True):
    """
    Given a solvated system (BSS object) and BSS protocol, run a process workflow with either
    Sander (CPU) or pmemd.cuda (GPU). NPT is typically done with GPU to save computing time.
    Returns the processed system.
    """

    # Create the process passing a working directory.
    if engine == "AMBER":
        if not pmemd:
            process = BSS.Process.Amber(system, protocol)
        elif pmemd:
            if AMBER_path is None:
                raise ValueError(
                    "AMBER path must be specified to run using AMBER engine."
                )
            process = BSS.Process.Amber(
                system,
                protocol,
                exe=AMBER_path,
                is_gpu=True,
            )
    elif engine == "GROMACS":
        process = BSS.Process.Gromacs(system, protocol)
    elif engine == "OpenMM":
        process = BSS.Process.OpenMM(system, protocol)
    # Start the process.
    process.start()

    # Wait for the process to exit.
    process.wait()

    # Check for errors.
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise _Exceptions.ThirdPartyError("The process exited with an error!")

    # If it worked, try to get the system. No need to block, since it's already finished.
    system = process.getSystem()

    return system


# In[ ]:


# Initialise the node object
node = BSS.Gateway.Node("Equilibrate a solvated ligand.")
# Set the node author and license.
node.addAuthor(
    name="Julien Michel",
    email="julien.michel@ed.ac.uk",
    affiliation="University of Edinburgh",
)
node.addAuthor(
    name="Matthew Burman",
    email="matthew@openbiosim.org",
    affiliation="OpenBioSim",
)
node.setLicense("GPLv3")

### Set the node inputs
node.addInput(
    "file",
    BSS.Gateway.FileSet(
        help="A topology/coordinate representation of a solvated ligand."
    ),
)

node.addInput(
    "nvt_restrained",
    BSS.Gateway.Time(
        help="The duration of the NVT restrained equilibration stage.",
        default=5,
        unit="picoseconds",
    ),
)
node.addInput(
    "nvt",
    BSS.Gateway.Time(
        help="The duration of the NVT restrained equilibration stage.",
        default=50,
        unit="picoseconds",
    ),
)
node.addInput(
    "npt",
    BSS.Gateway.Time(
        help="The duration of the NVT restrained equilibration stage.",
        default=200,
        unit="picoseconds",
    ),
)
node.addInput(
    "MDengine",
    BSS.Gateway.String(
        help="The MD engine to use for equilibration/minimisation.",
        allowed=["AMBER", "GROMACS", "OpenMM"],
        default="AMBER",
    ),
)
node.addInput(
    "AMBER_path",
    BSS.Gateway.String(
        help="Path to AMBER installation",
        default=None,
    ),
)
node.addInput(
    "output suffix",
    BSS.Gateway.String("Suffix of the output file", default="_equilibrated"),
)

node.addInput(
    "output directory",
    BSS.Gateway.String(
        "Name of the directory in which to save the output",
        default="equilibrated_free_systems",
    ),
)

node.addInput(
    "file_prefix",
    BSS.Gateway.String(
        help="Prefix for output files. Required for any nodes run in multiple instances."
    ),
)
### Set the node outputs
node.addOutput("system_eq", BSS.Gateway.FileSet(help="The equilibrated system."))

node.showControls()


engine = node.getInput("MDengine")
print("The simulation engine is %s" % engine)

name = Path(node.getInput("file")[0]).stem.split(".")[0]
outpath = Path("./" + node.getInput("output directory"))
# make full name of output
output_name_full = str(outpath / (name + node.getInput("output suffix")))

# lets first check if the full output files already exist, if they do then there is nothing to do
if (
    Path(output_name_full + ".prm7").exists()
    and Path(output_name_full + ".rst7").exists()
):
    print("Output files already exist, node exiting")
    # if they do we will not run the node again
    node.setOutput(
        "minimised", [output_name_full + ".prm7", output_name_full + ".rst7"]
    )
    node.validate(file_prefix=node.getInput("file_prefix"))
    exit()


outpath.mkdir(parents=True, exist_ok=True)

#######################################
### Load the system  ##
#######################################
system = BSS.IO.readMolecules(node.getInput("file"))


print(f"NVT equilibration while restraining all non-solvent atoms..")
protocol = BSS.Protocol.Equilibration(
    runtime=node.getInput("nvt_restrained"),
    temperature_start=0 * BSS.Units.Temperature.kelvin,
    temperature_end=300 * BSS.Units.Temperature.kelvin,
    restraint="all",
)
equil1 = runProcess(
    system, protocol, engine=engine, AMBER_path=node.getInput("AMBER_path")
)


# In[ ]:


print(f"NVT equilibration  without restraints..")
protocol = BSS.Protocol.Equilibration(
    runtime=node.getInput("nvt"), temperature=300 * BSS.Units.Temperature.kelvin
)

equil2 = runProcess(equil1, protocol, engine=engine)


# In[ ]:


print(f"NPT equilibration while restraining non-solvent heavy atoms..")
protocol = BSS.Protocol.Equilibration(
    runtime=node.getInput("npt"),
    pressure=1 * BSS.Units.Pressure.atm,
    temperature=300 * BSS.Units.Temperature.kelvin,
    restraint="heavy",
)
equil3 = runProcess(
    equil2, protocol, engine=engine, AMBER_path=node.getInput("AMBER_path")
)


# In[ ]:


print(f"NPT equilibration  without restraints..")
protocol = BSS.Protocol.Equilibration(
    runtime=node.getInput("npt"),
    pressure=1 * BSS.Units.Pressure.atm,
    temperature=300 * BSS.Units.Temperature.kelvin,
)
system_eq = runProcess(
    equil3, protocol, engine=engine, AMBER_path=node.getInput("AMBER_path")
)


# In[ ]:


# Save systems
node.setOutput(
    "system_eq", BSS.IO.saveMolecules(output_name_full, system_eq, ["PRM7", "RST7"])
)


# In[ ]:


node.validate(file_prefix=node.getInput("file_prefix"))


# In[ ]:
