#!/bin/sh

GENOME=GRCh38
TEMPLATE=slurm_cellranger_aggr
SAMPLENAME=All_conditions

mkdir -p scripts

WORKDIR=$(pwd)

# build libraries.csv file
(echo "library_id,fragments,cells"
for LIBRARY in $(cat All_conditions.txt); do
  echo "${LIBRARY},${WORKDIR}/${LIBRARY}/outs/fragments.tsv.gz,${WORKDIR}/${LIBRARY}/outs/singlecell.csv"
done) >	${SAMPLENAME}.csv

sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXWORKDIRXX%$WORKDIR% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/aggr_${SAMPLENAME}.slurm
