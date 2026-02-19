#!/bin/bash

OUTDIR=$1

for file in $OUTDIR/snap*.hdf5; do
  sbatch /home/r408l055/run-gen-pk-file.sh "$file" $OUTDIR
done
