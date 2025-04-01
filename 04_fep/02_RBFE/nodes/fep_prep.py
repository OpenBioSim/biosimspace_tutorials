import BioSimSpace as BSS
from pathlib import Path

from BioSimSpace.Types import Time
from BioSimSpace import _Exceptions


def prepare_fep_free(
    ligand1_name: str,
    ligand2_name: str,
    ligand1_files: list,
    ligand2_files: list,
    output_location: str | Path,
    num_lambda: int = 11,
    lambda_values: list = None,
    md_engine: str = "GROMACS",
    runtime: str | Time = 2 * BSS.Units.Time.nanosecond,
    setup_only: bool = True,
) -> bool:

    if md_engine not in ["GROMACS", "SOMD", "SOMD2"]:
        raise ValueError(
            "The md_engine must be either 'GROMACS', 'SOMD' or 'SOMD2'. Please check your input."
        )

    if lambda_values:
        if len(lambda_values) != num_lambda:
            raise ValueError(
                "The number of lambda values provided does not match the number of lambda points."
            )

    # first lets make the output directory
    outpath = Path(output_location)
    outpath_full = outpath / f"{ligand1_name}~{ligand2_name}" / "free"
    outpath_full.mkdir(parents=True, exist_ok=True)
    # print(f"Working directory: {outpath_full}")

    # load ligand files to BSS
    ligand_1_sys = BSS.IO.readMolecules(ligand1_files)
    ligand_2_sys = BSS.IO.readMolecules(ligand2_files)

    ligand_1 = ligand_1_sys[0]
    ligand_2 = ligand_2_sys[0]

    # print("Mapping and aligning...")
    mapping = BSS.Align.matchAtoms(ligand_1, ligand_2, complete_rings_only=True)
    inv_mapping = {v: k for k, v in mapping.items()}
    ligand_2_a = BSS.Align.rmsdAlign(ligand_2, ligand_1, inv_mapping)

    # Generate merged molecule.
    # print("Merging..")
    merged_ligs = BSS.Align.merge(ligand_1, ligand_2_a, mapping)

    ligand_1_sys.removeMolecules(ligand_1)
    ligand_1_sys.addMolecules(merged_ligs)
    system_free = ligand_1_sys

    if isinstance(runtime, str):
        try:
            runtime = Time(runtime)
        except ValueError as e:
            raise ValueError(
                f"Could not convert runtime {runtime} to a Time object. Please provide a valid time string: {e}"
            )

    freenrg_protocol = BSS.Protocol.FreeEnergyProduction(
        num_lam=num_lambda, runtime=runtime, lam_vals=lambda_values
    )

    if md_engine == "SOMD2":
        import yaml as _yaml

        if not lambda_values:
            somd2_config = {
                "num_lambda": num_lambda,
                "runtime": str(runtime),
                "output_directory": str(outpath_full),
            }
        else:
            somd2_config = {
                "num_lambda": num_lambda,
                "runtime": str(runtime),
                "lambda_values": lambda_values,
                "output_directory": str(outpath_full),
            }
    # now we will need 2 separate sets of processes, one for setup only and one for a full run
    if setup_only:
        # start with SOMD
        if md_engine == "SOMD":
            _ = BSS.FreeEnergy.Relative(
                system_free,
                freenrg_protocol,
                engine="SOMD",
                work_dir=str(outpath_full),
                setup_only=True,
            )
            del _
        if md_engine == "SOMD2":
            _ = BSS.Stream.save(
                system_free,
                str(outpath_full / f"{ligand1_name}~{ligand2_name}"),
            )
            with open(outpath_full / "config.yaml", "w") as yaml_file:
                _yaml.dump(somd2_config, yaml_file)
            del _

        if md_engine == "GROMACS":
            _ = BSS.FreeEnergy.Relative(
                system_free,
                freenrg_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full),
                setup_only=True,
            )
            del _
    # Now the more complex case in which we actually want to run a production simulation
    else:
        # print(f"Running FEP using {md_engine}")
        if md_engine == "SOMD":
            process = BSS.FreeEnergy.Relative(
                system_free,
                freenrg_protocol,
                engine="SOMD",
                work_dir=str(outpath_full),
            )
            process.run()
            process.wait()

        if md_engine == "SOMD2":
            # print("SOMD2 not currently supported for non setup runs")
            pass
            # #print(type(system_free._sire_object))
            # import somd2

            # somd2_config = somd2.config.Config(**somd2_config)

            # runner = somd2.runner.Runner(system_free._sire_object, somd2_config)
            # runner.run()

        if md_engine == "GROMACS":
            # We will need some additional min/eq for GROMACS
            min_protocol = BSS.Protocol.FreeEnergyMinimisation(
                num_lam=num_lambda, lam_vals=lambda_values
            )
            nvt_protocol = BSS.Protocol.FreeEnergyEquilibration(
                num_lam=num_lambda, pressure=None, lam_vals=lambda_values
            )
            npt_protocol = BSS.Protocol.FreeEnergyEquilibration(
                num_lam=num_lambda,
                pressure=1 * BSS.Units.Pressure.atm,
                lam_vals=lambda_values,
            )

            # print("Running minimisation...")
            FE_min = BSS.FreeEnergy.Relative(
                system_free,
                min_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full / "min"),
            )
            FE_min.run()
            FE_min.wait()

            # print("Running heating...")
            FE_heat = BSS.FreeEnergy.Relative(
                system_free,
                nvt_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full / "heat"),
            )
            FE_heat.run()
            FE_heat.wait()

            # print("Running equilibration...")
            FE_eq = BSS.FreeEnergy.Relative(
                system_free,
                npt_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full / "eq"),
            )
            FE_eq.run()
            FE_eq.wait()

            # print("Running production...")
            FE_prod = BSS.FreeEnergy.Relative(
                system_free,
                freenrg_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full),
            )
            FE_prod.run()
            FE_prod.wait()
    return True


def prepare_fep_bound(
    ligand1_name: str,
    ligand2_name: str,
    ligand1_files: list,
    ligand2_files: list,
    output_location: str | Path,
    num_lambda: int = 11,
    lambda_values: list = None,
    md_engine: str = "GROMACS",
    runtime: str | Time = 2 * BSS.Units.Time.nanosecond,
    setup_only: bool = True,
) -> None:

    if md_engine not in ["GROMACS", "SOMD", "SOMD2"]:
        raise ValueError(
            "The md_engine must be either 'GROMACS', 'SOMD' or 'SOMD2'. Please check your input."
        )

    if lambda_values:
        if len(lambda_values) != num_lambda:
            raise ValueError(
                "The number of lambda values provided does not match the number of lambda points."
            )

    # first lets make the output directory
    outpath = Path(output_location)
    outpath_full = outpath / f"{ligand1_name}~{ligand2_name}" / "bound"
    outpath_full.mkdir(parents=True, exist_ok=True)
    # print(f"Working directory: {outpath_full}")

    # load ligand files to BSS
    system_1 = BSS.IO.readMolecules(ligand1_files)
    system_2 = BSS.IO.readMolecules(ligand2_files)

    # Extract ligands and protein. Do this based on nAtoms and nResidues, as sometimes
    # the order of molecules is switched, so we can't use index alone.
    # bugfix in BSS makes the below redundant but keeping this in to be 100% sure we're getting the correct structures.
    system_ligand_1 = None
    protein = None
    n_residues = [mol.nResidues() for mol in system_1]
    n_atoms = [mol.nAtoms() for mol in system_1]
    for i, (n_resi, n_at) in enumerate(zip(n_residues[:20], n_atoms[:20])):
        if n_resi == 1 and n_at > 5:
            system_ligand_1 = system_1.getMolecule(i)
        elif n_resi > 1:
            protein = system_1.getMolecule(i)
        else:
            pass

    # loop over molecules in system to extract the ligand
    system_ligand_2 = None

    n_residues = [mol.nResidues() for mol in system_2]
    n_atoms = [mol.nAtoms() for mol in system_2]
    for i, (n_resi, n_at) in enumerate(zip(n_residues, n_atoms)):
        # grab the system's ligand and the protein. ignore the waters.
        if n_resi == 1 and n_at > 5:
            system_ligand_2 = system_2.getMolecule(i)
        else:
            pass

    if system_ligand_1 and system_ligand_2 and protein:
        # print("Using molecules ligand_1, ligand_2, protein:")
        # print(system_ligand_1, system_ligand_2, protein)
        pass
    else:
        raise _Exceptions.AlignmentError(
            "Could not extract ligands or protein from input systems. Check that your ligands/proteins are properly prepared by BSSligprep.sh!"
        )

    # Align ligand2 on ligand1
    # print("Mapping..")
    mapping = BSS.Align.matchAtoms(
        system_ligand_1, system_ligand_2, complete_rings_only=True
    )
    inv_mapping = {v: k for k, v in mapping.items()}

    # print("Aligning..")
    system_ligand_2_a = BSS.Align.rmsdAlign(
        system_ligand_2, system_ligand_1, inv_mapping
    )

    # Generate merged molecule.
    # print("Merging..")
    system_merged_ligs = BSS.Align.merge(system_ligand_1, system_ligand_2_a, mapping)

    system_1.removeMolecules(system_ligand_1)
    system_1.addMolecules(system_merged_ligs)
    system_bound = system_1

    if isinstance(runtime, str):
        try:
            runtime = Time(runtime)
        except ValueError as e:
            raise ValueError(
                f"Could not convert runtime {runtime} to a Time object. Please provide a valid time string: {e}"
            )

    freenrg_protocol = BSS.Protocol.FreeEnergyProduction(
        num_lam=num_lambda, runtime=runtime, lam_vals=lambda_values
    )

    if md_engine == "SOMD2":
        import yaml as _yaml

        if not lambda_values:
            somd2_config = {
                "num_lambda": num_lambda,
                "runtime": str(runtime),
                "output_directory": str(outpath_full),
            }
        else:
            somd2_config = {
                "num_lambda": num_lambda,
                "runtime": str(runtime),
                "lambda_values": lambda_values,
                "output_directory": str(outpath_full),
            }
    # now we will need 2 separate sets of processes, one for setup only and one for a full run
    if setup_only:
        # start with SOMD
        if md_engine == "SOMD":
            BSS.FreeEnergy.Relative(
                system_bound,
                freenrg_protocol,
                engine="SOMD",
                work_dir=str(outpath_full),
                setup_only=True,
            )
        if md_engine == "SOMD2":
            BSS.Stream.save(
                system_bound,
                str(outpath_full / f"{ligand1_name}~{ligand2_name}"),
            )
            with open(outpath_full / "config.yaml", "w") as yaml_file:
                _yaml.dump(somd2_config, yaml_file)

        if md_engine == "GROMACS":
            BSS.FreeEnergy.Relative(
                system_bound,
                freenrg_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full),
                setup_only=True,
            )
    # Now the more complex case in which we actually want to run a production simulation
    else:
        # print(f"Running FEP using {md_engine}")
        if md_engine == "SOMD":
            process = BSS.FreeEnergy.Relative(
                system_bound,
                freenrg_protocol,
                engine="SOMD",
                work_dir=str(outpath_full),
            )
            process.run()
            process.wait()

        if md_engine == "SOMD2":
            # print("SOMD2 not currently supported for non setup runs")
            pass
            # #print(type(system_free._sire_object))
            # import somd2

            # somd2_config = somd2.config.Config(**somd2_config)

            # runner = somd2.runner.Runner(system_free._sire_object, somd2_config)
            # runner.run()

        if md_engine == "GROMACS":
            # We will need some additional min/eq for GROMACS
            min_protocol = BSS.Protocol.FreeEnergyMinimisation(
                num_lam=num_lambda, lam_vals=lambda_values
            )
            nvt_protocol = BSS.Protocol.FreeEnergyEquilibration(
                num_lam=num_lambda, pressure=None, lam_vals=lambda_values
            )
            npt_protocol = BSS.Protocol.FreeEnergyEquilibration(
                num_lam=num_lambda,
                pressure=1 * BSS.Units.Pressure.atm,
                lam_vals=lambda_values,
            )
            # print("Running minimisation...")
            FE_min = BSS.FreeEnergy.Relative(
                system_bound,
                min_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full / "min"),
            )
            FE_min.run()
            FE_min.wait()

            # print("Running heating...")
            FE_heat = BSS.FreeEnergy.Relative(
                system_bound,
                nvt_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full / "heat"),
            )
            FE_heat.run()
            FE_heat.wait()

            # print("Running equilibration...")
            FE_eq = BSS.FreeEnergy.Relative(
                system_bound,
                npt_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full / "eq"),
            )
            FE_eq.run()
            FE_eq.wait()

            # print("Running production...")
            FE_prod = BSS.FreeEnergy.Relative(
                system_bound,
                freenrg_protocol,
                engine="GROMACS",
                work_dir=str(outpath_full),
            )
            FE_prod.run()
            FE_prod.wait()
    return True


if __name__ == "__main__":
    # should not be run as main
    raise RuntimeError("This script is not intended to be run as a standalone script.")
