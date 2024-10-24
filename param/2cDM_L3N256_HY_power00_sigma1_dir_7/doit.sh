#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=8g
#SBATCH --time=6:00:00
##SBATCH --mail-user=rtlow@ku.edu
##SBATCH --mail-type=ALL
#SBATCH --job-name=2cDM_L3N256_HY_power00_sigma1_dir_7
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib

# Which code directory
i=7
# Is this a restart run?
RESTART=0
# Point to output directory
OUTDIR="/home/r408l055/scratch/output"
# Point to final snapshot
FINAL_SNAP="${OUTDIR}/2cDM_L3N256_HY_power00_sigma1_dir_7/snap_007.hdf5"


# Modification is rarely needed
RESTART_DIR="${OUTDIR}/2cDM_L3N256_HY_power00_sigma1_dir_7/restartfiles/"
JOBNAME=$SLURM_JOB_NAME
PARAM='2cDM_L3N256_HY_power00_sigma1_dir_7.txt'
PARAM_PATH="../RUNS/boxes/${PARAM}"
LOG_PATH=/home/r408l055/scratch/logs/LOG_${JOBNAME}_$(date +"%Y_%m_%d_%H_%M_%S")

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
