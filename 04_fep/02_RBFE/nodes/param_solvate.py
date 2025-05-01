# A node for parameterising and solvating a single system
import BioSimSpace as BSS

from pathlib import Path

from uuid import uuid4

node = BSS.Gateway.Node("Parameterise and solvate a single free leg system")


# Set the node author and license.
node.addAuthor(
    name="Matthew Burman", email="matthew@openbiosim.org", affiliation="OpenBioSim"
)
node.setLicense("GPLv3")

node.addInput("file", BSS.Gateway.File(help="full path to the ligand file"))

node.addInput(
    "protein files",
    BSS.Gateway.FileSet(help="full path to parameterised protein files (rst7 and prm7"),
)

node.addInput(
    "ligand forcefield",
    BSS.Gateway.String(
        help="The name of the force field to use for parameterisation.",
        allowed=BSS.Parameters.forceFields(),
        default="gaff2",
    ),
)

node.addInput(
    "water model",
    BSS.Gateway.String(
        help="The name of the water model to use for solvation.",
        allowed=BSS.Solvent.waterModels(),
        default="tip3p",
    ),
)

node.addInput(
    "box length",
    BSS.Gateway.Length(
        help="Length of the box edges. Applied in x, y and z directions",
        unit="nanometer",
        default=5,
    ),
)

node.addInput(
    "box type",
    BSS.Gateway.String(
        help="Box type to use for both bound and free legs",
        allowed=BSS.Box.boxTypes(),
        default="cubic",
    ),
)

node.addInput(
    "ion conc",
    BSS.Gateway.Float(
        help="The ionic concentration in mol/litre.", minimum=0, maximum=1, default=0
    ),
)

node.addInput(
    "output suffix",
    BSS.Gateway.String(help="Suffix of the output file", default="solvated"),
)

node.addInput(
    "output directory",
    BSS.Gateway.String(
        "Name of the directory in which to save the output",
        default="param_solv_systems",
    ),
)

node.addInput(
    "file_prefix",
    BSS.Gateway.String(
        help="Prefix for output files. Required for any nodes run in multiple instances.",
        default="output",
    ),
)

node.addOutput(
    "free solvated",
    BSS.Gateway.FileSet(
        help="Files containing the parameterised and solvated free leg system."
    ),
)
node.addOutput(
    "bound solvated",
    BSS.Gateway.FileSet(
        help="Files containing the parameterised and solvated bound leg system."
    ),
)

node.showControls()

name = Path(node.getInput("file")).stem.split(".")[0]
outpath = Path("./" + node.getInput("output directory"))
# make full name of output
output_name_full_free = str(
    outpath / (name + "_free_" + node.getInput("output suffix"))
)
output_name_full_bound = str(
    outpath / (name + "_bound_" + node.getInput("output suffix"))
)

# lets first check if the full output files already exist, if they do then there is nothing to do
if (
    Path(output_name_full_free + ".prm7").exists()
    and Path(output_name_full_free + ".rst7").exists()
    and Path(output_name_full_bound + ".prm7").exists()
    and Path(output_name_full_bound + ".rst7").exists()
):
    print("Output files already exist, node exiting")
    # if they do we will not run the node again
    node.setOutput(
        "free solvated",
        [output_name_full_free + ".prm7", output_name_full_free + ".rst7"],
    )
    node.setOutput(
        "bound solvated",
        [output_name_full_bound + ".prm7", output_name_full_bound + ".rst7"],
    )
    node.validate(file_prefix=node.getInput("file_prefix"))
    exit()


system = BSS.IO.readMolecules(node.getInput("file"))

outpath.mkdir(parents=True, exist_ok=True)

# parameterise the ligand using the input forcefield
paramed = BSS.Parameters.parameterise(system[0], node.getInput("ligand forcefield"))
lig_p = paramed.getMolecule()


protein_folder = Path("../inputs/")
# protein is already parameterised
protein = BSS.IO.readMolecules(node.getInput("protein files"))
system = lig_p + protein

box_min, box_max = lig_p.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min, box_max)]
box_sizes = [x + node.getInput("box length") for x in box_size]

boxtype_dict = {
    "cubic": BSS.Box.cubic,
    "rhombicDodecahedronHexagon": BSS.Box.rhombicDodecahedronHexagon,
    "rhombicDodecahedronSquare": BSS.Box.rhombicDodecahedronSquare,
    "truncatedOctahedron": BSS.Box.truncatedOctahedron,
}
boxtype_func = boxtype_dict[node.getInput("box type")]


box, angles = boxtype_func(max(box_sizes))
# print(f"box of free leg: {box}, angles: {angles}")
lig_p_solvated = BSS.Solvent.solvate(
    node.getInput("water model"), molecule=lig_p, box=box, angles=angles, ion_conc=0.15
)
# again for p-l complex
box_min, box_max = system.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min, box_max)]
box_sizes = [x + node.getInput("box length") for x in box_size]

box, angles = boxtype_func(max(box_sizes))
# print(f"box of bound leg: {box}, angles: {angles}")
system_solvated = BSS.Solvent.solvate(
    node.getInput("water model"), molecule=system, box=box, angles=angles, ion_conc=0.15
)


node.setOutput(
    "free solvated",
    BSS.IO.saveMolecules(output_name_full_free, lig_p_solvated, ["PRM7", "RST7"]),
)

node.setOutput(
    "bound solvated",
    BSS.IO.saveMolecules(output_name_full_bound, system_solvated, ["PRM7", "RST7"]),
)

node.validate(file_prefix=node.getInput("file_prefix"))
