#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=ana
#SBATCH -o ../slurm_logs/ana_%A.out
#SBATCH -e ../slurm_logs/ana_%A.err

# sourcing
source $BSS
source $amber
source $gromacs
source $scripts_dir/extract_execution_model_bash.sh

# cd /home/anna/BioSimSpace
# git checkout feature-amber-pre-2023
# export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev"
# source $BSS
# cd $MAINDIRECTORY

date
echo "Folder for these runs is : $MAINDIRECTORY"
echo "Analysis for the transformation $1, $2."

python $scripts_dir/analysis_single_trans.py $1 $2
