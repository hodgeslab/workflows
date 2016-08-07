#!/bin/sh

WORKDIR=XXWORKDIRXX/regularizedFragLibrary
GENOME=mm9
FRAGMENTLENGTH=200

PATH=$PATH:~/src/bedtools2/bin
GENOMEFILE=~/data/UCSC_Downloads/${GENOME}/chrom.sizes

cd $WORKDIR

for BEDFILE in *.fragments_${FRAGMENTLENGTH}bp.bed; do
  PEAKBEDFILE=`echo $BEDFILE | sed s/.fragments_${FRAGMENTLENGTH}bp.bed/.peaks.bed/`
  OUTPUTFILE=`echo $BEDFILE | sed s/.bed/.bedGraph/`
  BGNAME=`echo $BEDFILE | sed s/.fragments_${FRAGMENTLENGTH}bp.bed//`
  BIGWIGFILE=`basename $OUTPUTFILE | sed s/.bedGraph/.bigWig/`

  if [ ! -r ${BIGWIGFILE} ]; then
    READCOUNT=`bedtools intersect -sorted -a $PEAKBEDFILE -b $BEDFILE | wc -l | awk '{ print 1e6/$1 }'`
    bedtools genomecov -scale $READCOUNT -bg -trackline -trackopts "name=\"$BGNAME\"" -i $BEDFILE -g $GENOMEFILE > $OUTPUTFILE
#    gunzip ${OUTPUTFILE}.gz
    bedGraphToBigWig $OUTPUTFILE $GENOMEFILE $BIGWIGFILE
    gzip $OUTPUTFILE
  fi
done

for i in `ls *.peaks.bed|cut -f 1 -d .|sed s/_rep.//|sort|uniq|tr "\n" " "`; do
  if [ ! -r ${i}_mean.fragments_${FRAGMENTLENGTH}bp.bigWig ]; then
    wiggletools mean ${i}_rep1.fragments_${FRAGMENTLENGTH}bp.bigWig ${i}_rep2.fragments_${FRAGMENTLENGTH}bp.bigWig > ${i}_mean.fragments_${FRAGMENTLENGTH}bp.wig
    wigToBigWig ${i}_mean.fragments_${FRAGMENTLENGTH}bp.wig $GENOMEFILE ${i}_mean.fragments_${FRAGMENTLENGTH}bp.bigWig
    rm ${i}_mean.fragments_${FRAGMENTLENGTH}bp.wig
  fi
done

