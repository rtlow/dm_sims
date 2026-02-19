#!/bin/bash

# Loop through the directories
for i in {0..9}
do
    dir=/home/r408l055/projects/SIDM/code$i 
    # Check if it is a directory
    if [ -d "$dir" ]; then
        # Enter the directory
        cd "$dir" || exit
        
        # submit job
        sbatch doit.sh

        # Return to the parent directory
        cd ..
    fi
done
