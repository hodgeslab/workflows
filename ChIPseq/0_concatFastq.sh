#!/bin/sh

FASTQDIR=./fastq
OUTDIR=./fastq
INFILE=fastqList.txt

mkdir -p $OUTDIR

while read -r line; do
  OUTFILE=$(echo "$line" | cut -f 1)
  PATTERN=$(echo "$line" | cut -f 2)
  echo "Assembling $OUTFILE"
  cat ${FASTQDIR}/${PATTERN} > ${OUTDIR}/${OUTFILE}
done < $INFILE
