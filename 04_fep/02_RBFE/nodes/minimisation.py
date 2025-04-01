import BioSimSpace as BSS
from pathlib import Path


def Minimisation(system, steps=10000, engine="AMBER"):
    protocol = BSS.Protocol.Minimisation(steps=steps)
    if engine == "GROMACS":
        process = BSS.Process.Gromacs(
            system, protocol, ignore_warnings=True, extra_args={"--ntmpi": 1}
        )
    elif engine == "AMBER":
        process = BSS.Process.Amber(
            system,
            protocol,
            is_gpu=True,
            exe="/home/matthew/AMBER/amber24/bin/pmemd.cuda",
        )
    else:
        raise TypeError("No valid MD engine")
    process.start()
    process.wait()
    # Check for errors.
    if process.isError():
        print(process.stdout())
        print(process.stderr())
    system = process.getSystem()
    return system


node = BSS.Gateway.Node(
    "A node used to minimise a single system. Assumes that inputs are in the form of a rst7 and parm7 file."
)


node.addAuthor(
    name="Matthew Burman",
    email="matthew@openbiosim.org",
    affiliation="OpenBioSim",
)
node.setLicense("GPLv3")


node.addInput("file", BSS.Gateway.FileSet(help="A rst7 and parm7 file pair."))
node.addInput(
    "steps",
    BSS.Gateway.Integer(
        help="The number of minimisation steps.",
        minimum=0,
        maximum=1000000,
        default=10000,
    ),
)

node.addInput(
    "output suffix",
    BSS.Gateway.String("Suffix of the output file", default="_minimised"),
)

node.addInput(
    "output directory",
    BSS.Gateway.String(
        "Name of the directory in which to save the output",
        default="minimised_systems",
    ),
)

node.addInput(
    "file_prefix",
    BSS.Gateway.String(
        help="Prefix for output files. Required for any nodes run in multiple instances."
    ),
)

node.addInput(
    "MDengine",
    BSS.Gateway.String(
        help="The MD engine to use for minimisation.",
        allowed=["AMBER", "GROMACS"],
        default="AMBER",
    ),
)

node.addOutput("minimised", BSS.Gateway.FileSet(help="The minimised molecular system"))

node.showControls()


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

system = BSS.IO.readMolecules(node.getInput("file"))

try:
    s_min = Minimisation(
        system, steps=node.getInput("steps"), engine=node.getInput("MDengine")
    )
except Exception as e:
    print("Error in minimisation")
    print(e)
    node.validate(file_prefix=node.getInput("file_prefix"))
    exit()

node.setOutput(
    "minimised",
    BSS.IO.saveMolecules(output_name_full, s_min, ["prm7", "rst7"]),
)

node.validate(file_prefix=node.getInput("file_prefix"))
