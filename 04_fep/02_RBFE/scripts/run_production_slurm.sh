#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --job-name=prod
#SBATCH -o ../slurm_logs/prod_%A_%a.out
#SBATCH -e ../slurm_logs/prod_%A_%a.err

# sourcing
module load cuda/11.6
module load gromacs/22.2
somdfreenrg="/export/users/XXX/sire.app/bin/somd-freenrg"
source $BSS
source extract_execution_model_bash.sh

# keep_traj="None"
# remove any trailing  '\r' characters
transformation=$(echo "$1" | tr -d '\r')
engine=$(echo "$2" | tr -d '\r')
num_lambda=$(echo "$3" | tr -d '\r')
date
echo "Folder for these runs is : $MAINDIRECTORY"
echo "The transformation is $transformation using $num_lambda windows and $engine as the MD engine"

# define no of windows based on sys arg 2
if [ $num_lambda = 11 ]; then
lamvals=( 0.0000 0.1000 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 )
fi
if [ $num_lambda = 17 ]; then
lamvals=( 0.0000 0.0625 0.1250 0.1875 0.2500 0.3125 0.3750 0.4375 0.5000 0.5625 0.6250 0.6875 0.7500 0.8125 0.8750 0.9375 1.0000 )
fi

# define lambda based on slurm array
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

# change to the trans dir, abort and message if not there
cd "$MAINDIRECTORY"/outputs/"$engine"/"$transformation" || exit
if [[ ! -d $MAINDIRECTORY/outputs/$engine/$transformation ]]; then
    echo "$MAINDIRECTORY/outputs/$engine/$transformation does not exist. Production run aborted..."
    exit 1
fi
trans_dir=$(pwd)

for dir in 'bound' 'free'; do

    if [ "$engine" = "SOMD" ]; then
        cd $dir
        cd lambda_$lam
        $somdfreenrg -c somd.rst7 -t somd.prm7 -m somd.pert -C somd.cfg -p CUDA

        cd $trans_dir

    elif [ "$engine" = "SOMD2" ]; then
        cd $dir
        somd2 "$transformation".bss --config config.yaml

        cd $trans_dir

    elif [ "$engine" = "GROMACS" ]; then

        cd $dir

        echo "min"
        gmx grompp -f min/lambda_$lam/gromacs.mdp -c min/lambda_$lam/initial_gromacs.gro -p min/lambda_$lam/gromacs.top -o min/lambda_$lam/gromacs.tpr
        gmx mdrun -ntmpi 1 -deffnm min/lambda_$lam/gromacs

        echo "heat"
        gmx grompp -f heat/lambda_$lam/gromacs.mdp -c min/lambda_$lam/gromacs.gro -p heat/lambda_$lam/gromacs.top -o heat/lambda_$lam/gromacs.tpr
        gmx mdrun -ntmpi 1 -deffnm heat/lambda_$lam/gromacs

        echo "eq"
        gmx grompp -f eq/lambda_$lam/gromacs.mdp -c heat/lambda_$lam/gromacs.gro -p eq/lambda_$lam/gromacs.top -t heat/lambda_$lam/gromacs.cpt -o eq/lambda_$lam/gromacs.tpr
        gmx mdrun -ntmpi 1 -deffnm eq/lambda_$lam/gromacs

        echo "prod"
        gmx grompp -f lambda_$lam/gromacs.mdp -c eq/lambda_$lam/gromacs.gro -p lambda_$lam/gromacs.top -t eq/lambda_$lam/gromacs.cpt -o lambda_$lam/gromacs.tpr
        gmx mdrun -ntmpi 1 -deffnm lambda_$lam/gromacs

    fi

    cd $trans_dir

done

echo "done."
