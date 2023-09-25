#!/usr/bin/env bash
echo $SLURM_JOB_NAME
date

lig0=$1
lig1=$2
engine=$3

sleep 5

echo "$BSSHOME/bin/python scripts/BSSanalyseFEP.py $lig0 $lig1 $engine"
$BSSHOME/bin/python scripts/BSSanalyseFEP.py $lig0 $lig1 $engine

exit 0
