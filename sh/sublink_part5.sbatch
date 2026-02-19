#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=6:00:00
#SBATCH --mail-user=rtlow@ku.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=sublink_part5
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib
#SBATCH -o /home/r408l055/scratch/logs/%x.%j.out

source ~/modules.sh
module load openmp

export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# Recommended memory:
# L75n1820FP, Galaxies: 30000 MB
# L75n910FP,  Subhalos: 7500 MB
# L75n1820FP, Subhalos: 60000 MB

TRACKING=Galaxies

BASEDIR=$1
TREEDIR=$1/sublink/$TRACKING
SNAPNUM_FIRST=0
SNAPNUM_LAST=127

# Create extended trees
echo ""; echo "CREATING EXTENDED TREES..."
python ~/sublink/Python/create_extended_trees.py $BASEDIR $TREEDIR $SNAPNUM_LAST

# Concatenate minimal trees
echo ""; echo "CONCATENATING MINIMAL TREES..."
python ~/sublink/Python/concatenate_trees.py ${TREEDIR}/tree

# Concatenate extended trees
echo ""; echo "CONCATENATING EXTENDED TREES..."
python ~/sublink/Python/concatenate_trees.py ${TREEDIR}/tree_extended

# Compute offsets
echo ""; echo "COMPUTING OFFSETS..."
python ~/sublink/Python/compute_offsets.py $BASEDIR $TREEDIR $SNAPNUM_FIRST $SNAPNUM_LAST

