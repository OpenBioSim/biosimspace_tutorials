#!/usr/bin/env bash
date

lig0=$1
lig1=$2
lambdastring=$3
engine=$4

IFS=',' read -r -a lambdas <<< "$lambdastring"

nwin=${#lambdas[@]}

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

echo "Array ID is" $ARRAY_TASK_ID


if [ "$ARRAY_TASK_ID" -ge "$nwin" ]; then
    stage="free"
    idx=$(( $ARRAY_TASK_ID - $nwin ))
else
    stage="bound"
    idx=$ARRAY_TASK_ID
fi


echo "idx is " $idx
lambda=${lambdas[$idx]}
echo "lambda is: "$lambda
echo "stage is: "$stage
echo "engine is: "$engine


if [[ $engine == *"SOMD"* ]]; then
    echo "USING SOMD"
    echo "cd outputs/SOMD/$lig0~$lig1/$stage/lambda_$lambda/"
    cd outputs/SOMD/$lig0~$lig1/$stage/lambda_$lambda/
    # run SOMD simulation.
    echo "$BSSHOME/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert 1> somd.log 2> somd.err"
    $BSSHOME/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert  1> somd.log 2> somd.err
    if [[ $REMOVE_SIM_DATA == "True" ]]; then
	    echo "Removing simulation data. If this is unintended, set REMOVE_SIM_DATA to False in processFEP-*.sh."
	    rm *.dcd *.s3* latest*
    fi

elif [[ $engine == *"GROMACS"* ]]; then
    echo "USING GROMACS"
    echo "cd outputs/GROMACS/$lig0~$lig1/$stage/lambda_$lambda/"
    cd outputs/GROMACS/$lig0~$lig1/$stage/lambda_$lambda/
    # run GROMACS simulation.
    echo "gmx mdrun -v -deffnm gromacs 1> gromacs.log 2> gromacs.err"
    gmx mdrun -v -deffnm gromacs 1> gromacs.log 2> gromacs.err
else
    echo "The FEP engine $engine is not supported"
fi

exit 0
