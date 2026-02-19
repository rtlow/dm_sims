#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: $O source_directory"
  exit 1
fi

source_dir="$1"

if [ ! -d "$source_dir" ]; then
  echo "Directory does not exist: $source_dir"
  exit 1
fi

for file in "$source_dir"/bak-*; do
  if [ -e "$file" ]; then
    filename=$(basename "$file")
    new_filename="${filename#bak-}"
    mv "$file" "$source_dir/$new_filename"
    echo "Moved: $file -> $source_dir/$new_filename"
  fi
done

echo "'bak-' restartfiles rewritten in $source_dir"
