#!/bin/sh

GENOME=mm10
TEMPLATE=slurm_cellranger_aggr
SAMPLENAME=All_conditions

mkdir -p scripts

WORKDIR=$(pwd)

# build libraries.csv file
(echo "library_id,molecule_h5"
for LIBRARY in $(cat All_conditions.txt); do
  echo "${LIBRARY},${WORKDIR}/${LIBRARY}/outs/molecule_info.h5"
done) >	${SAMPLENAME}.csv

sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXWORKDIRXX%$WORKDIR% $TEMPLATE > scripts/aggr_${SAMPLENAME}.slurm
