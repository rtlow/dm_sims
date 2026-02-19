#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=1g
#SBATCH --time=6:00:00
#SBATCH --mail-user=rtlow@ku.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=sublink_part2
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib
#SBATCH -o /home/r408l055/scratch/logs/%x.%j.out

source ~/modules.sh
module load openmpi

## OMP_NUM_THREADS controls the number of threads your application uses
## This variable cat be set by the following command:
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

TRACKING=Galaxies

# ----------------------- FIND DESCENDANTS

BASEDIR=$1

WRITEDIR=$1/sublink/$TRACKING/descendants
if [ ! -d $WRITEDIR ]; then
  mkdir -p $WRITEDIR
fi
WRITEPATH=$WRITEDIR/desc

SKIP_FILENAME=~/skip_snaps.txt

SNAPNUM_FIRST=0
SNAPNUM_LAST=127
SNAPNUM_START=0
SNAPNUM_END=127

# First pass
# PASS=first
# ~/sublink/Descendants/find_descendants $BASEDIR ${WRITEPATH}_${PASS} $SNAPNUM_FIRST $SNAPNUM_LAST $SNAPNUM_START $SNAPNUM_END $TRACKING ${PASS}

#~ # Second pass
PASS=second
~/sublink/Descendants/find_descendants $BASEDIR ${WRITEPATH}_${PASS} $SNAPNUM_FIRST $SNAPNUM_LAST $SNAPNUM_START $SNAPNUM_END $TRACKING ${PASS} ${SKIP_FILENAME}
#~ 
#~ # Compare descendants
#~ ~/sublink/Descendants/compare_descendants ${WRITEPATH} $SNAPNUM_FIRST $SNAPNUM_LAST ${SKIP_FILENAME}
#~ 
#~ # ---------------------------- BUILD TREES
#~ 
#~ TREEDIR=/n/ghernquist/vrodrigu/MergerTrees/output/$TRACKING/Illustris/${SIMULATION}_test
#~ INPUT_PATH=$TREEDIR/descendants/desc${SKIP_SUFFIX}
#~ OUTPUT_PATH=$TREEDIR/tree${SKIP_SUFFIX}
#~ 
#~ ~/sublink/SubhaloTrees/build_trees ${INPUT_PATH} ${OUTPUT_PATH} ${SNAPNUM_FIRST} ${SNAPNUM_LAST} ${SKIP_FILENAME}

sbatch sublink_part3.sbatch $1