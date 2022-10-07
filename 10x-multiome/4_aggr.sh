#!/bin/sh

GENOME=mm10
TEMPLATE=slurm_cellranger_aggr
SAMPLENAME=All_conditions

mkdir -p scripts

WORKDIR=$(pwd)

# build libraries.csv file
(echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info"
for LIBRARY in $(cat All_conditions.txt); do
  echo "${LIBRARY},${WORKDIR}/${LIBRARY}/outs/atac_fragments.tsv.gz,${WORKDIR}/${LIBRARY}/outs/per_barcode_metrics.csv,${WORKDIR}/${LIBRARY}/outs/gex_molecule_info.h5"
done) >	${SAMPLENAME}.csv

sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% \
    -e s%XXWORKDIRXX%$WORKDIR% \
    -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/aggr_${SAMPLENAME}.slurm
