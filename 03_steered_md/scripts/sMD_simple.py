import BioSimSpace as BSS
import os
from argparse import ArgumentParser
from shutil import copyfile


def parse_range(residue_range):
    """Parse residue range into a list

    Parameters
    ----------
    residue_range : str
        residue range, e.g. 179-184

    Returns
    -------
    residues : [int]
        a list of all residues in the range
    """
    start = int(residue_range.split("-")[0])
    end = int(residue_range.split("-")[1])

    residues = [i for i in range(start, end + 1)]

    return residues


def create_CV(reference, residues, system):
    """Create a collective variable to use in sMD

    Parameters
    ----------
    reference : BioSimSpace._SireWrappers._molecule.Molecule
        reference structure for RMSD calculation

    residues : [int]
        a list of residues that will be used to calculate RMSD

    system : BioSimSpace._SireWrappers._system.System
        steered system

    Returns
    -------
    rmsd_cv : BioSimSpace.Metadynamics.CollectiveVariable._rmsd.RMSD
        CV for use in steering
    """
    rmsd_indices = []
    for residue in reference.getResidues():
        if residue.index() in residues:
            for atom in residue.getAtoms():
                if atom.element() != "Hydrogen (H, 1)":
                    rmsd_indices.append(atom.index())

    rmsd_cv = BSS.Metadynamics.CollectiveVariable.RMSD(system, reference, rmsd_indices)

    return rmsd_cv


def setup_sMD_protocol(steering_runtime, total_runtime, force_constant, cv):
    """Create steered MD protocol

    Parameters
    ----------
    steering_runtime : BioSimSpace.Units.Time
        steering time. First 4 ps will be used to apply force

    total_runtime : BioSimSpace.Units.Time
        total simulation time. Usually steering_runtime + some short duration for relaxation of system

    force_constant : int
        force constant to use for steering

    Returns
    -------
    protocol : BioSimSpace.Protocol._steering.Steering
        steering protocol
    """
    start = 0 * BSS.Units.Time.nanosecond
    apply_force = 4 * BSS.Units.Time.picosecond

    nm = BSS.Units.Length.nanometer
    restraint_1 = BSS.Metadynamics.Restraint(cv.getInitialValue(), 0)
    restraint_2 = BSS.Metadynamics.Restraint(cv.getInitialValue(), force_constant)
    restraint_3 = BSS.Metadynamics.Restraint(0 * nm, force_constant)
    restraint_4 = BSS.Metadynamics.Restraint(0 * nm, 0)

    protocol = BSS.Protocol.Steering(
        cv,
        [start, apply_force, steering_runtime, total_runtime],
        [restraint_1, restraint_2, restraint_3, restraint_4],
        runtime=total_runtime,
        report_interval=2500,
        restart_interval=2500,
    )

    return protocol


def run_sMD(
    topology,
    coordinates,
    reference,
    residues,
    steering_runtime,
    total_runtime,
    force_constant,
):
    """Run steered MD with PLUMED and BioSimSpace.

    Parameters
    ----------
    topology : str
        AMBER topology file

    coordinates : str
        equilibrated system coordinates

    reference : str
        steering target PDB structure

    residues : [int]
        a list of residues that will be used to calculate RMSD

    steering_runtime : BioSimSpace.Units.Time
        steering time

    total_runtime : BioSimSpace.Units.Time
        total simulation runtime

    force_constant : int
        force to use for steering

    Returns
    -------
    None
    """
    system = BSS.IO.readMolecules([topology, coordinates])
    reference = BSS.IO.readMolecules(reference).getMolecule(0)

    rmsd_cv = create_CV(reference, residues, system)

    protocol = setup_sMD_protocol(
        steering_runtime, total_runtime, force_constant, rmsd_cv
    )

    process = BSS.Process.Amber(
        system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda'
    )

    process.start()
    process.wait()

    files = ["amber.nc", "COLVAR", "stdout", "amber.out"]
    for file in files:
        copyfile(f"{process.workDir()}/{file}", file)

    return None


def __main__():
    parser = ArgumentParser(
        description="Set up and run sMD with PLUMED and BioSimSpace"
    )
    parser.add_argument(
        "--topology", type=str, required=True, help="AMBER topology file"
    )
    parser.add_argument(
        "--coordinates", type=str, required=True, help="equilibrated system coordinates"
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="target molecule PDB"
    )
    parser.add_argument(
        "--residues",
        type=str,
        required=True,
        help="residues that will be used to calculate RMSD, given as a range, e.g. 179-185",
    )
    parser.add_argument(
        "--steering_runtime", type=float, required=True, help="steering time in ns"
    )
    parser.add_argument(
        "--total_runtime", type=float, required=True, help="total simulation runtime"
    )
    parser.add_argument(
        "--force", type=int, required=True, help="steering force_constant"
    )
    args = parser.parse_args()

    run_sMD(
        args.topology,
        args.coordinates,
        args.reference,
        parse_range(args.residues),
        args.steering_runtime * BSS.Units.Time.nanosecond,
        args.total_runtime * BSS.Units.Time.nanosecond,
        args.force,
    )

    return None


if __name__ == "__main__":
    __main__()
