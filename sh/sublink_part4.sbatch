#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=4g
#SBATCH --time=6:00:00
#SBATCH --mail-user=rtlow@ku.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=sublink_part4
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib
#SBATCH -o /home/r408l055/scratch/logs/%x.%j.out
#SBATCH --array=0-7

source ~/modules.sh
module load openmpi

export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# Recommended memory:
# L75n1820FP, Galaxies: 30000 MB
# L75n910FP,  Subhalos: 7500 MB
# L75n1820FP, Subhalos: 60000 MB

# ----------- General-purpose variables --------------

NUM_JOBS=8
TRACKING=Galaxies

BASEDIR=$1
TREEDIR=$1/sublink/$TRACKING
SNAPNUM_FIRST=0
SNAPNUM_LAST=127

# Create extra columns
echo ""; echo "CREATING EXTRA COLUMNS..."
python ~/sublink/Python/create_columns.py $BASEDIR $TREEDIR $SNAPNUM_FIRST $SNAPNUM_LAST $NUM_JOBS ${SLURM_ARRAY_TASK_ID}

sbatch sublink_part5.sbatch $1