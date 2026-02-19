#!/bin/bash

# Loop through the directories
for i in {0..9}
do
    dir=/home/r408l055/projects/SIDM/code$i
    # Check if it is a directory
    if [ -d "$dir" ]; then
        
        sbatch /home/r408l055/compile-dir.sh $dir

        echo Submitted job to compile $dir

    fi
done
