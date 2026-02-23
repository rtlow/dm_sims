#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=1g
#SBATCH --time=6:00:00
#SBATCH --mail-user=rtlow@ku.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=rockstar-L3N256
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib
#SBATCH -o /home/r408l055/scratch/logs/%x.%j.out

source ~/modules.sh

module load conda

conda activate cosmo

pip install mpi4py --user

INDIR=$1

mpiexec -np $SLURM_NTASKS python ~/rockstar-galaxies/rockstar-submit-L3N256.py $INDIR

sbatch ~/run-consistent-dir.sh $INDIR/Rockstar
