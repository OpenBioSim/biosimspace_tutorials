#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --job-name=ligprep
#SBATCH -o ../slurm_logs/ligprep_%A_%a.out
#SBATCH -e ../slurm_logs/ligprep_%A_%a.err

# sourcing
source $BSS
source $amber
source $gromacs
source extract_execution_model_bash.sh

date
echo "Folder for these runs is : $MAINDIRECTORY"
echo "Ligands file is : $lig_file"

lig=${lig_array[$SLURM_ARRAY_TASK_ID]}

echo "prep for $lig..."
python $scripts_dir/ligprep.py $lig
