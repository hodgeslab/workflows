#!/bin/sh

GENOME=mm10
EXPERIMENTS="BM_vehicle BM_BRM014"

module load hdf5
module load R

for EXPERIMENT in ${EXPERIMENTS}; do
  cd $EXPERIMENT

  # DoubletFinder
  if [ ! -r qc_DoubletFinder.R ]; then
    cp ../qc_DoubletFinder.R .
  fi

  if [ ! -r qc_DoubletFinder_pass_singlets.txt ]; then
    Rscript qc_DoubletFinder.R
  fi

  cd ..
done
