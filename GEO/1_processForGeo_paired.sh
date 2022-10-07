#!/bin/sh

INFILE=$1

FASTQFOLDER=../fastq
COUNTSFOLDER=""
BIGWIGFOLDER=../fragLibrary
FRAGMENTLEN=200

cat $INFILE | sed '$a\' | while read -r line; do
  PUBLICNAME=$(echo "$line" | cut -f 1)
  LOCALNAME=$(echo "$line" | cut -f 2)
  FASTQFILE1=$(echo "$line" | cut -f 3)
  FASTQFILE2=$(echo "$line" | cut -f 4)

  ln -s ${FASTQFOLDER}/$FASTQFILE1 ${PUBLICNAME}_R1.fastq.gz
  ln -s ${FASTQFOLDER}/$FASTQFILE2 ${PUBLICNAME}_R2.fastq.gz

  if [ ! -z $COUNTSFOLDER ] && [ -r ${COUNTSFOLDER}/${LOCALNAME}_counts.txt ]; then
    ln -s ${COUNTSFOLDER}/${LOCALNAME}_counts.txt ${PUBLICNAME}_counts.txt
  fi

  if [ ! -z $BIGWIGFOLDER ] && [ -r ${BIGWIGFOLDER}/${LOCALNAME}.fragments_${FRAGMENTLEN}bp.bigWig ]; then
    ln -s ${BIGWIGFOLDER}/${LOCALNAME}.fragments_${FRAGMENTLEN}bp.bigWig ${PUBLICNAME}.bigWig
  fi

done
