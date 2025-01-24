#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=8g
#SBATCH --time=6:00:00
##SBATCH --mail-user=rtlow@ku.edu
##SBATCH --mail-type=ALL
#SBATCH --job-name=2cDM_L3N512_DM_power00_sigma1_dir_5
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib

# Which code directory
i=5
# Is this a restart run?
RESTART=0
# Point to output directory
OUTDIR="/home/r408l055/scratch/output"

# Modification is rarely needed
JOBNAME=$SLURM_JOB_NAME
PARAM='2cDM_L3N512_DM_power00_sigma1_dir_5.txt'
PARAM_PATH="../RUNS/boxes/${PARAM}"
LOG_PATH=/home/r408l055/scratch/logs/LOG_${JOBNAME}_$(date +"%Y_%m_%d_%H_%M_%S")

# Actual run starting code
~/remove_core.sh

module load hdf5/1.10.5
module load openmpi/4.0

for j in $(seq 0 19);
do
  mpiexec ./Arepo $PARAM_PATH 3 $j >> $LOG_PATH
done

~/remove_core.sh
