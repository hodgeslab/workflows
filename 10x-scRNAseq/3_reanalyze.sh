#!/bin/sh

GENOME=mm10
TEMPLATE=slurm_cellranger_reanalyze
mkdir -p scripts

WORKDIR=$(pwd)
FASTQ=${WORKDIR}/fastq

cat sampleList.txt | sed '$a\' | while read -r line; do
  SAMPLENAME=$(echo "$line" | cut -f 1)
  LIBRARYID=$(echo "$line" | cut -f 2)
  LIBRARYTYPE=$(echo "$line" | cut -f 3)

  sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% \
      -e s%XXWORKDIRXX%$WORKDIR% \
      -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/reanalyze_${SAMPLENAME}.slurm

done
