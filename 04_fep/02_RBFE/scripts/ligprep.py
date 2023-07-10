import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import sys
import os
import csv

### Settings.
minim_steps = 250
runtime_short_nvt = 5  # ps
runtime_nvt = 50  # ps
runtime_npt = 200  # ps

main_dir = os.environ["MAINDIRECTORY"]

amber_home = os.environ["AMBERHOME"]
pmemd_path = amber_home + "/bin/pmemd.cuda"

protein_file = os.environ["protein_file"]
ligands_folder = os.environ["ligands_folder"]

################
### Open 'ligands.dat', find sys.argv[1] 'th entry
print(f"{sys.argv[0]} {sys.argv[1]}")
lig_name = sys.argv[1]

# Exit if prep files are already available for ligand and protein

if os.path.exists(f"{main_dir}/prep/protein/{lig_name}_sys_equil_solv.prm7"):
    if os.path.exists(f"{main_dir}/prep/ligands/{lig_name}_lig_equil_solv.prm7"):
        print(f" Prep files already generated for {lig_name}")
        sys.exit(0)
#################
### Load matching input with BSS.read.IO.
lig = BSS.IO.readMolecules(f"{ligands_folder}/{lig_name}.sdf")[0]

#################
### parameterise ligand by deriving requested FF from protocol.dat.
stream = open(f"{main_dir}/protocol.dat", "r")
lines = stream.readlines()
ff_query = lines[0].rstrip()

# do tests on protocol, parameterise with correct FF.
if not "ligand" in ff_query:
    raise NameError(
        "Please supply a ligand force field on the first line of protocol.dat. The line should look like (e.g.):\n"
        + "ligand forcefield = GAFF2"
    )
else:
    if "GAFF1" in ff_query:
        print(f"Parameterising {lig_name} using GAFF1 force field.")
        lig_p = BSS.Parameters.gaff(lig).getMolecule()

    elif "GAFF2" in ff_query:
        print(f"Parameterising {lig_name} using GAFF2 force field.")
        lig_p = BSS.Parameters.gaff2(lig).getMolecule()
    else:
        raise NameError(
            f"Force field not supported: {ff_query}. Please use either of [GAFF1, GAFF2, OpenForceField], or"
            + " edit this script to use other force fields available in BSS.Parameters.forceFields()."
        )
# at this point BSS should have thrown an error if parameterisation failed, so no checks needed.
### make a copy of ligand, solvate original ligand object.
lig_p_copy = lig_p.copy()

# need to derive some settings from protocol.dat again.
stream = open(f"{main_dir}/protocol.dat", "r")
lines = stream.readlines()

# get the solvent force field.
solvent_query = lines[2].rstrip().replace(" ", "").split("=")[-1]

# get the box size settings.
boxsize_query = lines[3].rstrip().replace(" ", "").split("=")[-1]
box_axis_length = boxsize_query.split("*")[0]
box_axis_unit = boxsize_query.split("*")[1]
if box_axis_unit.lower() == "nm" or box_axis_unit.lower() == "nanometer":
    box_axis_unit = BSS.Units.Length.nanometer

elif box_axis_unit.lower() == "a" or box_axis_unit.lower() == "angstrom":
    box_axis_unit = BSS.Units.Length.angstrom
else:
    raise NameError(
        "Input unit not recognised. Please use any of ['spc', 'spce', 'tip3p', 'tip4p', 'tip5p'] in "
        + "the fourth line of protocol.dat in the shape of (e.g.):\nbox edges = 10*angstrom"
    )


"""
figure out the ligand dimensions, add the specified water box together with the largest ligand axis.
This is a workaround to account for adding ions in some cases. Based on:
https://github.com/michellab/BioSimSpace/issues/111
"""
box_min, box_max = lig_p.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min, box_max)]
box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]


# do the same for ligand +protein system.
protein = BSS.IO.readMolecules([f"{protein_file}.rst7", f"{protein_file}.prm7"])[0]
system = lig_p + protein

# box size
boxsize_query = lines[3].rstrip().replace(" ", "").split("=")[-1]
box_axis_length = boxsize_query.split("*")[0]
box_axis_unit = boxsize_query.split("*")[1]
if box_axis_unit.lower() == "nm" or box_axis_unit.lower() == "nanometer":
    box_axis_unit = BSS.Units.Length.nanometer
elif box_axis_unit.lower() == "a" or box_axis_unit.lower() == "angstrom":
    box_axis_unit = BSS.Units.Length.angstrom
else:
    raise NameError(
        "Input unit not recognised. Please use any of ['spc', 'spce', 'tip3p', 'tip4p', 'tip5p'] in "
        + "the fourth line of protocol.dat in the shape of (e.g.):\nbox edges = 10*angstrom"
    )

# get the box type settings.
boxtype_query = lines[4].rstrip().replace(" ", "").split("=")[-1]

boxtype_dict = {
    "cubic": BSS.Box.cubic,
    "truncatedoctahedron": BSS.Box.truncatedOctahedron,
    "octahedral": BSS.Box.truncatedOctahedron,
}

if boxtype_query not in boxtype_dict:
    raise NameError(
        "Input box type not recognised. Please use any of ['orthorhombic', 'octahedral', 'triclinic']"
        + "in the fifth line of protocol.dat in the shape of (e.g.):\nbox type = orthorhombic"
    )


boxtype_func = boxtype_dict[boxtype_query]
print(f"Solvating for {lig_name}...")

box_min, box_max = lig_p.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min, box_max)]
box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]

box, angles = boxtype_func(max(box_sizes))
lig_p_solvated = BSS.Solvent.solvate(
    solvent_query, molecule=lig_p, box=box, angles=angles, ion_conc=0.15
)

box_min, box_max = system.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min, box_max)]
box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]

box, angles = boxtype_func(max(box_sizes))
system_solvated = BSS.Solvent.solvate(
    solvent_query, molecule=system, box=box, angles=angles, ion_conc=0.15
)

#################
### Minimise and equilibrate our systems to make the ligand relax inside the pocket.

# first save and reload the systems to bypass some amber issues with BSS waters.
BSS.IO.saveMolecules(
    f"{main_dir}/tmp/{lig_name}_lig_s", lig_p_solvated, ["PRM7", "RST7"]
)
BSS.IO.saveMolecules(
    f"{main_dir}/tmp/{lig_name}_sys_s", system_solvated, ["PRM7", "RST7"]
)

lig_p_solvated = BSS.IO.readMolecules(
    [f"{main_dir}/tmp/{lig_name}_lig_s.prm7", f"{main_dir}/tmp/{lig_name}_lig_s.rst7"]
)
system_solvated = BSS.IO.readMolecules(
    [f"{main_dir}/tmp/{lig_name}_sys_s.prm7", f"{main_dir}/tmp/{lig_name}_sys_s.rst7"]
)


def runProcess(system, protocol, engine="AMBER", pmemd=False):
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
            process = BSS.Process.Amber(system, protocol, exe=pmemd_path)

    elif engine == "GROMACS":
        process = BSS.Process.Gromacs(system, protocol)

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


############# first minimise/equilibrate the solvated ligand.
print("\n#### Working on solvated ligand.")

print(f"Minimising in {minim_steps} steps..")
protocol = BSS.Protocol.Minimisation(steps=minim_steps)
minimised = runProcess(lig_p_solvated, protocol)

print(
    f"PMEMD NVT equilibration for {runtime_short_nvt} ps while restraining all non-solvent atoms.."
)
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_short_nvt * BSS.Units.Time.picosecond,
    temperature_start=0 * BSS.Units.Temperature.kelvin,
    temperature_end=300 * BSS.Units.Temperature.kelvin,
    restraint="all",
)
equil1 = runProcess(minimised, protocol, pmemd=True)

print(f"PMEMD NVT equilibration for {runtime_nvt} ps without restraints..")
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_nvt * BSS.Units.Time.picosecond,
    temperature=300 * BSS.Units.Temperature.kelvin,
)

equil2 = runProcess(equil1, protocol, pmemd=True)

print(
    f"PMEMD NPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms.."
)
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_npt * BSS.Units.Time.picosecond,
    pressure=1 * BSS.Units.Pressure.atm,
    temperature=300 * BSS.Units.Temperature.kelvin,
    restraint="heavy",
)
equil3 = runProcess(equil2, protocol, pmemd=True)

print(f"PMEMD NPT equilibration for {runtime_npt} ps without restraints..")
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_npt * BSS.Units.Time.picosecond,
    pressure=1 * BSS.Units.Pressure.atm,
    temperature=300 * BSS.Units.Temperature.kelvin,
)
lig_equil_fin = runProcess(equil3, protocol, pmemd=True)


############# repeat for ligand + protein system. Include extra restrain="backbone" step.
print("\n#### Working on solvated ligand+protein.")
print(f"Minimising in {minim_steps} steps..")
protocol = BSS.Protocol.Minimisation(steps=minim_steps)
minimised = runProcess(system_solvated, protocol)


print(
    f"PMEMD NVT equilibration for {runtime_short_nvt} ps while restraining all non-solvent atoms.."
)
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_short_nvt * BSS.Units.Time.picosecond,
    temperature_start=0 * BSS.Units.Temperature.kelvin,
    temperature_end=300 * BSS.Units.Temperature.kelvin,
    restraint="all",
)
equil1 = runProcess(minimised, protocol, pmemd=True)

print(
    f"PMEMD NVT equilibration for {runtime_nvt} ps while restraining all backbone atoms.."
)
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_nvt * BSS.Units.Time.picosecond,
    temperature=300 * BSS.Units.Temperature.kelvin,
    restraint="backbone",
)
equil2 = runProcess(equil1, protocol, pmemd=True)

print(f"PMEMD NVT equilibration for {runtime_nvt} ps without restraints..")
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_nvt * BSS.Units.Time.picosecond,
    temperature_end=300 * BSS.Units.Temperature.kelvin,
)

equil3 = runProcess(equil2, protocol, pmemd=True)

print(
    f"PMEMD NPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms.."
)
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_npt * BSS.Units.Time.picosecond,
    pressure=1 * BSS.Units.Pressure.atm,
    temperature=300 * BSS.Units.Temperature.kelvin,
    restraint="heavy",
)
equil4 = runProcess(equil3, protocol, pmemd=True)

print(f"PMEMD NPT equilibration for {runtime_npt} ps without restraints..")
protocol = BSS.Protocol.Equilibration(
    runtime=runtime_npt * BSS.Units.Time.picosecond,
    pressure=1 * BSS.Units.Pressure.atm,
    temperature=300 * BSS.Units.Temperature.kelvin,
)
sys_equil_fin = runProcess(equil4, protocol, pmemd=True)

# finally, save last snapshot of both equilibrated objects.
os.system("mkdir -p prep/ligands")
os.system("mkdir -p prep/protein")

print("Saving solvated/equilibrated systems.")
print("\n Ligand:")
print(lig_equil_fin)
BSS.IO.saveMolecules(
    f"{main_dir}/prep/ligands/{lig_name}_lig_equil_solv",
    lig_equil_fin,
    ["PRM7", "RST7"],
)

print("\n Ligand + protein:")
print(sys_equil_fin)
BSS.IO.saveMolecules(
    f"{main_dir}/prep/protein/{lig_name}_sys_equil_solv",
    sys_equil_fin,
    ["PRM7", "RST7"],
)
print("Done.")
