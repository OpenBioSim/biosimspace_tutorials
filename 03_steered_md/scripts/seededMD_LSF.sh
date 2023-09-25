#!/bin/bash
#
#BSUB -J "seededMD[1-3]"
#BSUB -q "phase12_slough"
#BSUB -n 1
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -R "rusage[mem=1] span[ptile=8] select[type==any]"
#BSUB -M 4000000
#BSUB -cwd /work/slough_md12_scratch/ahardie/seededMD
#BSUB -outdir /work/slough_md12_scratch/ahardie/seededMD
#BSUB -o lsf.%J.out
#BSUB -e lsf.%J.err

mkdir snapshots/snapshot_$LSB_JOBINDEX
source /home/model/MD-SOFTWARE/BSSenv.bashrc

python 02_run_seededMD.py --snapshot snapshot_${LSB_JOBINDEX}.pdb --phosphate_params /home/e628835/Documents/BioSimSpaceTutorials/02_steered_md/data/phosphate_params/ --output snapshot_${LSB_JOBINDEX}/
