#!/bin/sh

GENOME=mm10
TEMPLATE=slurm_cellranger_reanalyze_aggr
mkdir -p scripts

WORKDIR=$(pwd)
SAMPLENAME=All_conditions
OUTDIR=${SAMPLENAME}_singlets
mkdir -p ${OUTDIR}

INDEX=1
for LIBRARY in $(cat All_conditions.txt); do
  cat ${LIBRARY}/qc_combined_pass_singlets.txt | sed "s/-1$/-${INDEX}/"
  INDEX=$(expr $INDEX + 1)
done > ${SAMPLENAME}/qc_combined_pass_singlets.txt

sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% \
    -e s%XXWORKDIRXX%$WORKDIR% \
    -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/reanalyze_aggr_${SAMPLENAME}.slurm
