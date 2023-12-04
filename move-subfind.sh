#!/bin/bash

# Destination directory where files will be copied
destination_directory="./fof-subhalo/"

# Create destination directory if it doesn't exist
mkdir -p "$destination_directory"

# Loop through each source directory provided as a command-line argument
for source_directory in "$@"; do
    # Find all files in subfolders starting with "my-starting-string" and copy them
    find "$source_directory" -type f -name 'fof_subhalo_tab*' -exec cp --parents {} "$destination_directory" \;
done
