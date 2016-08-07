#!/bin/sh

WORKDIR=XXWORKDIRXX
GENOME=XXGENOMEXX

PATH=$PATH:~/src/bedtools2/bin
GENOMEFILE=~/data/UCSC_Downloads/${GENOME}/chrom.sizes

cd $WORKDIR

for NAME in `cut -f 1 sampleList.txt`; do
  BAMFILE=${NAME}.bam
  BEDGRAPHFILE=${NAME}.bedGraph
  BIGWIGFILE=${NAME}.bigWig

  if [ ! -r ${BEDGRAPHFILE}.gz ]; then
    READCOUNT=`samtools view -c $BAMFILE | awk '{ print 1e6/$1 }'`
    bedtools genomecov -scale $READCOUNT -bg -split -trackline -trackopts "name=\"$NAME\"" -ibam $BAMFILE -g $GENOMEFILE > $BEDGRAPHFILE
#    gunzip ${BEDGRAPHFILE}.gz
    bedGraphToBigWig $BEDGRAPHFILE $GENOMEFILE $BIGWIGFILE
    gzip $BEDGRAPHFILE
  fi
done

for i in `cut -f 1 sampleList.txt|sed s/_rep.//|sort|uniq|tr "\n" " "`; do
  if [ ! -r ${i}_mean.bigWig ]; then
    wiggletools mean ${i}_rep1.bigWig ${i}_rep2.bigWig > ${i}_mean.wig
    wigToBigWig ${i}_mean.wig $GENOMEFILE ${i}_mean.bigWig
    rm ${i}_mean.wig
  fi
done

