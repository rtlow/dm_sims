#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8g
#SBATCH --time=6:00:00
#SBATCH --mail-user=rtlow@ku.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=compile
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib
#SBATCH -o /home/r408l055/scratch/logs/%x.%j.out
#SBATCH -e /home/r408l055/scratch/logs/%x.%j.err

source ~/modules.sh

cd "$1"

make clean && make -j$(nproc)
