#!/bin/bash

# Loop through the directories
for dir in /home/r408l055/scratch/output/*
do
    # Check if it is a directory
    if [ -d "$dir" ]; then
        if [ -f "$dir/snap_127.hdf5" ]; then
          sbatch /home/r408l055/run-rockstar-dir.sh $dir
        else
          echo "Final snapshot doesn't exist for $dir"
        fi
    fi
done
