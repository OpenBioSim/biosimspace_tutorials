#!/usr/bin/env bash


# retrieve the array ID depending on whether we're running on LSF or slurm.
if ((${#LSB_JOBINDEX[@]})); then
   echo "This is an LSF job."
   ARRAY_TASK_ID=$(($LSB_JOBINDEX-1))

elif ((${#SLURM_ARRAY_TASK_ID[@]})); then
   echo "This is a Slurm job."
   ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID

else
   echo "No array ID found - are you sure you are running this script on a cluster?"
fi

echo "Array ID is"
date
$BSSHOME/bin/ipython scripts/BSSligprep.py $ARRAY_TASK_ID

date

exit 0
