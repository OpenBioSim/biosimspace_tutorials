#!/bin/bash
#SBATCH --job-name=steeredMD
#SBATCH --gres=gpu:1
#SBATCH --time=80:00:00
#SBATCH -o steered_MD.out

#load modules
#ensure AMBERHOME is set
source /home/adele/software/amber22/amber.sh

#get the right python environment
source /home/adele/anaconda3/etc/profile.d/conda.sh
conda activate bss

#run python script
#output will be saved in current directory
python scripts/sMD_simple.py --topology data/system.prm7 --coordinates data/system.rst7 --reference data/reference.pdb --residues 178-184 --steering_runtime 150 --total_runtime 152 --force 3500

