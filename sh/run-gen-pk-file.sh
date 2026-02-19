#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=32g
#SBATCH --time=6:00:00
##SBATCH --mail-user=rtlow@ku.edu
##SBATCH --mail-type=ALL
#SBATCH --job-name=Gen-PK
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib
#SBATCH -o /home/r408l055/scratch/logs/%x.$j.out
#SBATCH -e /home/r408l055/scratch/logs/%x.%j.err

module load hdf5
module load fftw3

file=$1
OUTDIR=$2

/home/r408l055/GenPK/gen-pk -i "$file" -o $OUTDIR
