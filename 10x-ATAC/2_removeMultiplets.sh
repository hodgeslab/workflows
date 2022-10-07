#!/bin/sh

GENOME=hg38
EXPERIMENTS="A549_EV A549_SMARCA4"

module load hdf5
module load R

for EXPERIMENT in ${EXPERIMENTS}; do
  cd $EXPERIMENT

  # AMULET
  WORKDIR=$(pwd)
  AMULETHOME=/s1/opt/AMULET
  FRAGMENTS=outs/atac_fragments.tsv.gz
  BARCODES=outs/per_barcode_metrics.csv
  AUTOSOMES=${AMULETHOME}/${GENOME}_autosomes.txt
  REPEATS=/s1/share/UCSC_Downloads/${GENOME}/${GENOME}-blacklist.bed
  OUTPUTDIR=./

  cat $BARCODES | sed s/,is_cell,/,is__cell_barcode,/ > qc_AMULET_barcodes.txt
  BARCODES=qc_AMULET_barcodes.txt

  ${AMULETHOME}/AMULET.sh \
    $FRAGMENTS \
    $BARCODES \
    $AUTOSOMES \
    $REPEATS \
    $OUTPUTDIR \
    ${AMULETHOME}/

  cd ..
done
