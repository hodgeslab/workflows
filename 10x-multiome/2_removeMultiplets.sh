#!/bin/sh

GENOME=mm10
EXPERIMENTS="PB_vehicle BM_vehicle PB_BRM014 BM_BRM014"

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

  # combine and list only those passing both using R
  Rscript -e "singlets <- read.table('qc_DoubletFinder_pass_singlets.txt', header=F, stringsAsFactors=F)[,1];
              multiplets <- read.table('MultipletBarcodes_01.txt', header=F, stringsAsFactors=F)[,1];
              write.table(setdiff(singlets,multiplets), 'qc_combined_pass_singlets.txt', quote=F, row.names=F, col.names=F)"

  cd ..
done
