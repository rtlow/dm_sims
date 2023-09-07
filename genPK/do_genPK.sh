#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: $O <folder>"
  exit 1
fi

folder="$1"

if [ ! -d "$folder" ]; then
  echo "Folder $folder does not exist."
  exit 1
fi

outpath="/home/ryan/genPK/"
fext=".hdf5"

for i in {0..7}; do
  fname="snap_00${i}"
  fname+=$fext
  file="$folder/$fname"
  if [ -f "$file" ]; then
    ~/code/GenPK/gen-pk -i "$file" -o "$outpath"
  else
    echo "$file not found."
  fi
done
