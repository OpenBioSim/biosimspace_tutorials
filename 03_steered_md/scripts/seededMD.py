import BioSimSpace as BSS
import os
from shutil import copyfile
from argparse import ArgumentParser


def run_MD(system, runtime=50 * BSS.Units.Time.nanosecond, output_dir=""):
    """Run equilibrium MD with AMBER (pmemd.cuda). Requires $AMBERHOME to be set.

    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        prepared system

    runtime : BioSimSpace.Units.Time
        simulation duration

    output_dir : str
        output directory

    Returns
    -------
    None
    """
    protocol = BSS.Protocol.Production(runtime=runtime)
    process = BSS.Process.Amber(
        system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda'
    )

    process.start()
    process.wait()

    # copy over results
    [
        copyfile(f"{process.workDir()}/{file}", f"{output_dir}{file}")
        for file in ["amber.nc", "stdout", "amber.rst7"]
    ]

    return None


def seededMD(snapshot, topology, output_dir=""):
    """Run an equilibrium MD simulation using AMBER (pmemd.cuda). Requires $AMBERHOME to be set to a correct AMBER installation.

    Parameters
    ----------
    snapshot : str
        snapshot PDB path

    phosphate_parameters : str
        path to additional phospho residue parameters

    output_dir : str
        output directory

    Returns
    -------
    None
    """
    system = BSS.IO.readMolecules([snapshot, topology])
    run_MD(system, output_dir=output_dir)

    return None


def __main__():
    parser = ArgumentParser(
        description="Set up and run seeded MD from snapshots. For PTP1B with peptide substrate"
    )
    parser.add_argument(
        "--snapshot", type=str, required=True, help="path to snapshot coordinates"
    )
    parser.add_argument(
        "--topology", type=str, required=True, help="path to system topology"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="",
        help="output directory. If not specified, output will be saved in current directory",
    )
    args = parser.parse_args()

    seededMD(args.snapshot, args.topology, args.output)

    return None


if __name__ == "__main__":
    __main__()
