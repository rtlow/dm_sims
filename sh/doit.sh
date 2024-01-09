#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=8g
#SBATCH --time=6:00:00
##SBATCH --mail-user=rtlow@ku.edu
##SBATCH --mail-type=ALL
#SBATCH --job-name=2cDM_L3N256_power00_sigma10
#SBATCH --partition=sixhour
#SBATCH --constraint "intel"
#SBATCH --constraint=ib

# Which code directory
i=5
# Is this a restart run?
RESTART=0
# Point to output directory
OUTDIR="/home/r408l055/scratch/output"
# Point to final snapshot
FINAL_SNAP="${OUTDIR}/live${i}/snap_007.hdf5"


# Modification is rarely needed
RESTART_DIR="${OUTDIR}/live${i}/restartfiles/"
JOBNAME=$SLURM_JOB_NAME
PARAM="param_L3N256_${i}.txt"
PARAM_PATH="../RUNS/710/${PARAM}"
LOG_PATH='/home/r408l055/scratch/logs/LOG_${JOBNAME}_$(date +"%Y_%m_%d_%H_%M_%S")'

# Flags for iteration go/no-go
RESTART_EXISTS=false
FINAL_SNAP_EXISTS=false
RESTART_WRITTEN_SUCCESS=true

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
	FINAL_SNAP_EXISTS=true
	echo "Final snapshot for Job $JOBNAME was written. Bye!"
fi

# Check for the existence of restartfiles
for i in $(seq 0 $((N_RESTART - 1))); do
    rfile="restart.${i}"
	file_name="${RESTART_DIR}/restart.${i}"
    if [ ! -e "$file_name" ]; then
        echo "File $rfile does not exist."
    fi
	else
	    RESTART_EXISTS=true
done

# Check last half of logfile
# for successful restart write
total_lines=$(wc -l < $LOG_PATH)
nlines=$((total_lines / 2 + 1))
tail_output=$(tail -n $nlines $LOG_PATH)
search_strings=('RESTART: Writing restart files' 'RESTART: Backup restart files' 'RESTART: load/save took ' 'RESTART: done.')

for string in "${search_strings[@]}"; do
  if ! grep -q "$string" <<< "$tail_output"; then
    RESTART_WRITTEN_SUCCESS=false
	echo "Failed to rewrite restart files!"
	echo "Check for stale restart files!"
    break
  fi
done

# Iteration?
if [[ $RESTART_EXISTS && $RESTART_WRITTEN_SUCCESS && ! $FINAL_SNAP_EXISTS  ]]; then
  sbatch doit-restart.sh
else
  echo "Job $JOBNAME stopped."
fi
