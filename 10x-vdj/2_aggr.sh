#!/bin/sh

GENOME=mm10
TEMPLATE=slurm_cellranger_aggr
SAMPLENAME=All_conditions

mkdir -p scripts

WORKDIR=$(pwd)

# build libraries.csv file
(echo "sample_id,vdj_contig_info,donor,origin"
for LIBRARY in $(cat All_conditions.txt); do
  echo "${LIBRARY},${WORKDIR}/${LIBRARY}/outs/vdj_contig_info.pb,${LIBRARY},${LIBRARY}"
done) >	${SAMPLENAME}.csv

sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXWORKDIRXX%$WORKDIR% $TEMPLATE > scripts/aggr_${SAMPLENAME}.slurm
