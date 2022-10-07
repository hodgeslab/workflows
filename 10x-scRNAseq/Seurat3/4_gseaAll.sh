#!/bin/sh

WORKDIR=$(pwd)
TEMPLATE=slurm_gsea
EXPERIMENTS="1-2"

mkdir -p scripts

for EXPERIMENT in $EXPERIMENTS; do
  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXEXPERIMENTXX%$EXPERIMENT% $TEMPLATE > scripts/gsea_${EXPERIMENT}.slurm
done
