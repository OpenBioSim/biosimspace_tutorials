#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --job-name=fepprep
#SBATCH -o ../slurm_logs/fepprep_%A_%a.out
#SBATCH -e ../slurm_logs/fepprep_%A_%a.err

# sourcing
module load cuda/11.6
module load amber/22
module load gromacs/20.4
source $BSS
source extract_execution_model_bash.sh

date
echo "Folder for these runs is : $MAINDIRECTORY"
echo "Network file is : $net_file"

trans=${trans_array[$SLURM_ARRAY_TASK_ID]}
eng=${eng_array[$SLURM_ARRAY_TASK_ID]}
win=${win_array[$SLURM_ARRAY_TASK_ID]}

echo "fepprep for $trans using $eng"
python $scripts_dir/fepprep.py $trans
