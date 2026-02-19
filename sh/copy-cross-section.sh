#!/bin/bash

CROSS_FOLDER=/home/r408l055/projects/SIDM/cross_sections/$1

# Loop through the directories
for i in {0..9}
do
    dir=/home/r408l055/projects/SIDM/code$i
    # Check if it is a directory
    if [ -d "$dir" ]; then
        # Enter the directory
        cd "$dir" || exit
        
        rm sidm_cross_reaction*

        if [ -d "$CROSS_FOLDER" ]; then
          # Copy desired config files
          ln -s $CROSS_FOLDER/* .
          
          echo Linked $1 to $dir
        fi

        # Return to the parent directory
        cd ..
    fi
done
