#!/bin/bash

# Loop through the directories
for dir in /home/r408l055/scratch/to_process/*
do
    # Check if it is a directory
    if [ -d "$dir" ]; then
      sbatch /home/r408l055/run-rockstar-L3N256.sh $dir
    fi
done
