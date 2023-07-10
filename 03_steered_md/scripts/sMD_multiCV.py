from importlib.metadata import files
import BioSimSpace as BSS
import os
from argparse import ArgumentParser
from shutil import copyfile


def create_rmsd_cv(system, reference):
    """Create WPD loop RMSD CV for sMD.

    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        steered system

    reference : BioSimSpace._SireWrappers._molecule.Molecule
        reference structure for RMSD calculation

    Returns
    -------
    rmsd_cv : BioSimSpace.Metadynamics.CollectiveVariable._rmsd.RMSD
        CV for use in steering
    """
    rmsd_indices = []
    for residue in reference.getResidues():
        if 178 >= residue.index() <= 184:
            for atom in residue.getAtoms():
                if atom.element() != "Hydrogen (H, 1)":
                    rmsd_indices.append(atom.index())

    rmsd_cv = BSS.Metadynamics.CollectiveVariable.RMSD(system, reference, rmsd_indices)

    return rmsd_cv


def create_torsion_cv(system):
    """Create a CV for the Chi1 angle of Tyr152

    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        steered system

    Returns
    -------
    torsion_cv : BioSimSpace.Metadynamics.CollectiveVariable._torsion.Torsion
        CV for use in steering
    """
    torsion_indices = []
    for atom in system.getMolecule(0).getResidues()[152].getAtoms():
        if atom.name() in ["N", "CA", "CB", "CG"]:
            torsion_indices.append(atom.index())

    torsion_cv = BSS.Metadynamics.CollectiveVariable.Torsion(torsion_indices)

    return torsion_cv


def create_distance_cv(system):
    """Create a CV for the stacking of Phe196 to Phe280, defined as the distance
    between the CG atoms.

    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        steered system

    Returns
    -------
    distance_cv : BioSimSpace.Metadynamics.CollectiveVariable._distnace.Distance
        CV for use in steering
    """
    distance_indices = []
    for residue in system.getMolecule(0).getResidues():
        if residue.index() == 196 or residue.index() == 280:
            for atom in residue.getAtoms():
                if atom.name() == "CG":
                    distance_indices.append(atom.index())
                    break

    distance_cv = BSS.Metadynamics.CollectiveVariable.Distance(
        distance_indices[0], distance_indices[1]
    )

    return distance_cv


def create_restraints(force):
    """Create restraints for each CV at each steering step

    Parameters
    ----------
    force : int
        force constant

    Returns
    -------
    restraints : [[BioSimSpace.Metadynamics.Restraints]]
        steering restraints
    """
    nm = BSS.Units.Length.nanometer
    rad = BSS.Units.Angle.radian

    restraints = [
        [
            BSS.Metadynamics.Restraint(0.64 * nm, 0),
            BSS.Metadynamics.Restraint(1.1 * rad, 0),
            BSS.Metadynamics.Restraint(0.56 * nm, 0),
        ],  # initial
        [
            BSS.Metadynamics.Restraint(0.64 * nm, force),
            BSS.Metadynamics.Restraint(1.1 * rad, force),
            BSS.Metadynamics.Restraint(0.56 * nm, force),
        ],  # apply force
        [
            BSS.Metadynamics.Restraint(0 * nm, force),
            BSS.Metadynamics.Restraint(1.1 * rad, force),
            BSS.Metadynamics.Restraint(0.4 * nm, force),
        ],  # steering
        [
            BSS.Metadynamics.Restraint(0 * nm, 0),
            BSS.Metadynamics.Restraint(1.1 * rad, 0),
            BSS.Metadynamics.Restraint(0.4 * nm, 0),
        ],
    ]  # release force

    return restraints


def create_protocol(system, reference, steering_runtime, total_runtime, force):
    """Create steering protocol

    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        steered system

    reference : BioSimSpace._SireWrappers._molecule.Molecule
        reference structure for RMSD calculation

    steering_runtime : BioSimSpace.Units.Time
        steering time. First 4 ps will be used to apply force

    total_runtime : BioSimSpace.Units.Time
        total simulation time. Usually steering_runtime + some short duration for relaxation of system

    Returns
    -------
    protocol : BioSimSpace.Protocol._steering.Steering
        steering protocol
    """
    rmsd_cv = create_rmsd_cv(system, reference)
    torsion_cv = create_torsion_cv(system)
    distance_cv = create_distance_cv(system)

    restraints = create_restraints(force)

    ns = BSS.Units.Time.nanosecond
    schedule = [0 * ns, 0.04 * ns, steering_runtime, total_runtime]

    protocol = BSS.Protocol.Steering(
        [rmsd_cv, torsion_cv, distance_cv],
        schedule,
        restraints,
        runtime=total_runtime,
        report_interval=2500,
        restart_interval=2500,
    )

    return protocol


def run_sMD(topology, coordinates, reference, steering_runtime, total_runtime, force):
    """Run multi-CV steered MD, steering PTP1B from inactive conformation (open WPD loop) to
    active conformation (closed WPD loop). The CVs are predefined, but the force constant and
    steering duration can be changed.

    Parameters
    ----------
    topology : str
        topology file

    coordinates : str
        coordinate file

    reference : str
        RMSD reference PDB file

    steering_runtime : BioSimSpace.Units.Time
        steering duration

    total_runtime : BioSimSpace.Units.Time
        total simulation duration

    force : int
        force constant
    """
    system = BSS.IO.readMolecules([topology, coordinates])
    reference = BSS.IO.readMolecules(reference).getMolecule(0)

    protocol = create_protocol(
        system, reference, steering_runtime, total_runtime, force
    )

    process = BSS.Process.Amber(
        system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda'
    )
    process.start()
    process.wait()

    files = [
        "amber.nc",
        "COLVAR",
        "amber.rst7",
        "plumed.dat",
        "reference_1.pdb",
        "amber.out",
    ]

    for file in files:
        copyfile(f"{process.workDir()}/{file}", f"./{file}")

    return None


def __main__():
    parser = ArgumentParser(
        description="Run multi-CV steered MD, steering PTP1B from inactive conformation (open WPD loop) to active conformation (closed WPD loop)."
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
        args.steering_runtime * BSS.Units.Time.nanosecond,
        args.total_runtime * BSS.Units.Time.nanosecond,
        args.force,
    )

    return None


if __name__ == "__main__":
    __main__()
