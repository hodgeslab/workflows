#!/bin/sh

# Each line in sampleList.txt contains a single library with the tab-delimited format:
# samplename	bigWig	bedfile
# In this file, each replicate is on its own line. If multiple bedfiles are needed for
# each sample, then use multiple lines that each contain a separate bedfile.

WORKDIR=$(pwd)
TEMPLATE=slurm_matrix

mkdir -p scripts

cat sampleList.txt | sed '$a\' | while read -r line; do
  SAMPLENAME=$(echo "$line" | cut -f 1)
  BIGWIG=$(echo "$line" | cut -f 2)
  BEDFILE=$(echo "$line" | cut -f 3)
  BEDNAME=$(basename $BEDFILE .bed)

  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXSAMPLENAMEXX%${SAMPLENAME}% -e s%XXBIGWIGXX%$BIGWIG% -e s%XXBEDFILEXX%$BEDFILE% $TEMPLATE > scripts/matrix_${BEDNAME}_${SAMPLENAME}.slurm
done
