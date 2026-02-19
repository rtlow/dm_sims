#!/bin/bash

POWER=$1

SIGMAS=(
  0.5
  0.646
  0.834
  1.08
  1.39
  1.8
  2.32
  3.0
  3.87
  5.0
)

INDEX=0

CROSS_FOLDER=/home/r408l055/projects/SIDM/cross_sections/$1

CODE_FOLDERS=/home/r408l055/projects/SIDM/code*
# Loop through the directories
for dir in $CODE_FOLDERS
do
    # Check if it is a directory
    if [ -d "$dir" ]; then
        # Enter the directory
        cd "$dir" || exit
         
        rm sidm_cross_reaction*

        SIGMA="${SIGMAS[$INDEX]}"

        CROSS_FOLDER=/home/r408l055/projects/SIDM/cross_sections/${POWER}_sigma${SIGMA}
        
        if [ -d "$CROSS_FOLDER" ]; then

          # Copy desired config files
          ln -s $CROSS_FOLDER/* .

          echo Linked $CROSS_FOLDER to "$dir"
        fi
        
        let INDEX=${INDEX}+1
        # Return to the parent directory
        cd ..
    fi
done
