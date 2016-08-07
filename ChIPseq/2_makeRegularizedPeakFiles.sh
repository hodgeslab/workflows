#!/bin/sh

PATH=$PATH:~/src/bedtools2/bin
INDIR=.
OUTDIR=regularizedFragLibrary
FRAGMENTLENGTH=200
mkdir -p $OUTDIR

cut -f 1 sampleList.txt > All_conditions.txt

for i in `cat All_conditions.txt`; do
  FRAGMENTSBED=${i}.fragments_${FRAGMENTLENGTH}bp.bed
  cp ${INDIR}/$FRAGMENTSBED $OUTDIR

  BED1=${OUTDIR}/${i}.peaks.bed.tmp
  BED2=${OUTDIR}/${i}.peaks.bed
  bedtools intersect -sorted -c -a ${INDIR}/${i}_peaks.bed -b ${INDIR}/${i}.fragments_${FRAGMENTLENGTH}bp.bed > $BED1
  total=`awk 'BEGIN{FS="\t"; a=0}{a=a+$NF}END{print a/1e6}' $BED1`
  awk -v total=$total 'BEGIN{FS="\t"; OFS="\t"}{$5=$6;$6=$6/total;print}' $BED1 > $BED2
  rm ${BED1}
done

