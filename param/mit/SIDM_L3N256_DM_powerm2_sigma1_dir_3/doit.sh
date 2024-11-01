#!/bin/bash

# Job Name and Files
#SBATCH -J SIDM_L3N256_DM_powerm2_sigma1_dir_3

#Output and error
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory
#SBATCH -D ./

#Partition & time limit
#SBATCH --partition=sched_mit_mvogelsb
#SBATCH --time=120:00:00

#Number of nodes and MPI tasks per node:
#SBATCH --nodes=2

#SBATCH --mem=0 # get all memory on each node
#SBATCH --ntasks=128 # less tasks than cores available, hopefully will give more mem per cpu

#SBATCH --constraint=centos7

#SBATCH --exclusive

#SBATCH --export=ALL
#SBATCH --exclude=node1259

#SBATCH --mail-user=rtlow@ku.edu
#SBATCH --mail-type=ALL


# Which code directory
i=3
# Is this a restart run?
RESTART=0
# Point to output directory
OUTDIR="/home/rtlow/scratch/output"
# Point to final snapshot
FINAL_SNAP="${OUTDIR}/SIDM_L3N256_DM_powerm2_sigma1_dir_3/snap_007.hdf5"


# Modification is rarely needed
RESTART_DIR="${OUTDIR}/SIDM_L3N256_DM_powerm2_sigma1_dir_3/restartfiles/"
JOBNAME=$SLURM_JOB_NAME
PARAM='SIDM_L3N256_DM_powerm2_sigma1_dir_3.txt'
PARAM_PATH="../RUNS/boxes/${PARAM}"
LOG_PATH=/home/rtlow/scratch/logs/LOG_${JOBNAME}_$(date +"%Y_%m_%d_%H_%M_%S")

# Flags for iteration go/no-go
RESTART_EXISTS=1
N_FINAL_SNAP_EXISTS=0
N_RUN_TERMINATED=0
N_EXISTING_RUNS=0

GO=0
# Actual run starting code
~/remove_core.sh

module load hdf5/1.10.5
module load openmpi/4.0
mpiexec ./Arepo $PARAM_PATH $RESTART > $LOG_PATH

~/remove_core.sh

# Iteration go/no-go logic
N_RESTART=$(($SLURM_JOB_NUM_NODES*$SLURM_NTASKS_PER_NODE))

# Check if the final snapshot was written
if [ -f "$FINAL_SNAP" ]; then
	N_FINAL_SNAP_EXISTS=0
	echo "Final snapshot for Job $JOBNAME was written. Bye!"
else
  N_FINAL_SNAP_EXISTS=1
  echo "Still snapshots to go..."
fi

# Check for the existence of restartfiles
for i in $(seq 0 $((N_RESTART - 1))); do
  rfile="restart.${i}"
	file_name="${RESTART_DIR}/restart.${i}"
    if [ ! -e "$file_name" ]; then
        echo "File $rfile does not exist."
        RESTART_EXISTS=0
    fi
done

TERMINATE_STRING="TERMINATE: ******!!!!!******"

# Check if this run was successful
if grep -q "$TERMINATE_STRING" "$LOG_PATH"; then
  N_RUN_TERMINATED=0
  echo "Job $JOBNAME was TERMINATED!!!!"
  echo "Look for unintended output!"
else
  echo "Previous job good..."
  N_RUN_TERMINATED=1
fi

# Check if other jobs with this name are running
count=$( squeue --format="%30j" | grep -c "$JOBNAME" )

if [[ $count > 1 ]]; then
  N_EXISTING_RUNS=0
  echo "Other runs exist! NO go!!!"
else
  N_EXISTING_RUNS=1
  echo "No other runs."
fi
GO=$(($RESTART_EXISTS * $N_FINAL_SNAP_EXISTS * $N_RUN_TERMINATED * $N_EXISTING_RUNS ))

echo "Go for restart? $GO"

if (($GO)) ; then
  echo "GO GO GO"
  sbatch doit-restart.sh
else
  echo "Job $JOBNAME stopped."
fi
