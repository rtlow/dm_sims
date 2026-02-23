#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=6:00:00
#SBATCH --mail-user=rtlow@ku.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=consistent-trees
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib
#SBATCH -o /home/r408l055/scratch/logs/%x.%j.out

source /home/r408l055/modules.sh

INDIR=$1

perl ~/rockstar-galaxies/scripts/gen_merger_cfg.pl $INDIR/rockstar.cfg

cd ~/consistent-trees

perl do_merger_tree.pl $INDIR/outputs/merger_tree.cfg
