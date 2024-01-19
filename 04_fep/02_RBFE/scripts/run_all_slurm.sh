#! /bin/bash
#
# Run all the lig prep, FEP prep, and the production runs

# TODO set the main directory filepath, make sure BSS branch and engines are sourced correctly, replace all XXX
export MAINDIRECTORY="XXX"                             # Set file path ; this should be the absolute file path to the main folder
export scripts_dir="$MAINDIRECTORY/scripts"            # choose location of scripts
export protein_file="$MAINDIRECTORY/inputs/prot_water" # this should be the prm7 and rst7 file name (without the extension) for the protein
export ligands_folder="$MAINDIRECTORY/inputs/ligands"  # location of the input files for the ligands

module load cuda/11.6
module load amber/22
module load gromacs/22.2
export BSS="/export/users/XXX/anaconda3/bin/activate biosimspace-dev"
source $BSS

# might have to source like this, depending on where it is
# export amber="/usr/local/amber22/amber.sh"
# export gromacs="/usr/local/gromacs/bin/GMXRC"
# source $amber
# source $gromacs

# export all files for later scripts
export lig_file="$MAINDIRECTORY/ligands.dat"
export net_file="$MAINDIRECTORY/network.dat"
export prot_file="$MAINDIRECTORY/protocol.dat"

# sourcing - as needed in the othe sh scripts
source extract_execution_model_bash.sh

###########
echo "The folder for all these runs is $MAINDIRECTORY"
echo ${lig_array[@]}
echo ${trans_array[@]}
echo ${eng_array[@]}
echo ${win_array[@]}

# make output dir for slurm out and err files
if [[ ! -d ../slurm_logs ]]; then
	mkdir ../slurm_logs
fi

# make directory for tmp from lig prep
if [[ ! -d ../tmp ]]; then
	mkdir ../tmp
fi

# chmod all files so can be executed by sbatch.
chmod u+x run_ligprep_slurm.sh
chmod u+x run_fepprep_slurm.sh
chmod u+x run_production_slurm.sh

# Run the runs
# ligand prep
jidlig=$(sbatch --parsable --array=0-$((${#lig_array[@]} - 1)) run_ligprep_slurm.sh)
echo "ligand prep jobid is $jidlig"

# FEP prep
jidfep=$(sbatch --dependency=afterany:${jidlig} --parsable --array=0-$((${#trans_array[@]} - 1)) run_fepprep_slurm.sh)
echo "FEP prep jobid is $jidfep"

# Production runs for the transformation
for i in "${!trans_array[@]}"; do

	# if using SOMD2 do not submit an array job
	if [[ " ${eng_array[*]} " =~ "SOMD2" ]]; then
		jidprod=$(sbatch --dependency=afterany:"${jidfep}" --parsable run_production_slurm.sh "${trans_array[i]}" "${eng_array[i]}")
		echo "Production jobid for ${trans_array[i]}, ${eng_array[i]} is $jidprod"
	else
		jidprod=$(sbatch --dependency=afterany:"${jidfep}" --parsable --array=0-$(("${win_array[i]}" - 1)) run_production_slurm.sh "${trans_array[i]}" "${eng_array[i]}" ${win_array[i]} $repeats)
		echo "Production jobid for ${trans_array[i]}, ${eng_array[i]} is $jidprod"
	fi
	jidana=$(sbatch --dependency=afterany:${jidprod} --parsable $scripts_dir/run_analysis_slurm.sh ${trans_array[i]} ${eng_array[i]})
	echo "Analysis jobid for ${trans_array[i]} is $jidana"
done
