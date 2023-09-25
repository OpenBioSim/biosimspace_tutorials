#!/bin/bash
#
#BSUB -J sMD
#BSUB -q "phase12_slough"
#BSUB -n 1
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -R "rusage[mem=1] span[ptile=8] select[type==any]"
#BSUB -M 4000000
#BSUB -cwd /work/slough_md12_scratch/ahardie
#BSUB -outdir /work/slough_md12_scratch/ahardie
#BSUB -o lsf.%J.out
#BSUB -e lsf.%J.err

cd /work/slough_md12_scratch/ahardie
#source /home/model/MD-SOFTWARE/amber18-gnu-cu10/amber.sh
source /home/model/MD-SOFTWARE/plumed-2.6.1/sourceme.sh
source /home/model/MD-SOFTWARE/BSSenv-new.bashrc

/home/e628835/software/sire.app/bin/python 01_run_sMD.py --topology system.prm7 --coordinates system.rst7 --reference reference.pdb --residues 179-185 --steering_runtime 150 --total_runtime 152 --force 3500
