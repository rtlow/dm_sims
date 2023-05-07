#!/bin/bash

# Job Name and Files
#SBATCH -J 2cDM_L3N512_power00_sigma1

#Output and error
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory
#SBATCH -D ./

#Partition & time limit
#SBATCH --partition=sched_mit_mvogelsb
#SBATCH --time=120:00:00

#Number of nodes and MPI tasks per node:
#SBATCH --nodes=8

#SBATCH --mem=0 # get all memory on each node
#SBATCH --ntasks=256 # less tasks than cores available, hopefully will give more mem per cpu

#SBATCH --constraint=centos7

#SBATCH --exclusive

#SBATCH --export=ALL
#SBATCH --exclude=node1259

#SBATCH --mail-user=rtlow@ku.edu
#SBATCH --mail-type=ALL



OUTDIR="/home/rtlow/output"
JOBNAME="2cDM_L3N512_DM_power00_sigma1"
FILE="${OUTDIR}/live2/snap_004.hdf5"

source /home/rtlow/modules.sh

mpiexec --mca opal_warn_on_missing_libcuda 0 \
       --mca btl '^openib' --mca pml ucx \
       -np $SLURM_NTASKS ./Arepo ../RUNS/710/param_L3N512_2.txt 0 > /home/rtlow/scratch/logs/LOG_${JOBNAME}_$(date +"%Y_%m_%d_%H_%M_%S")
       
#if [ ! -f "$FILE" ]; then
#  sbatch doit-restart.sh
#fi
