#!/bin/sh

WORKDIR=XXWORKDIRXX
GENOME=XXGENOMEXX

GENOMEFILE=/s1/share/UCSC_Downloads/${GENOME}/chrom.sizes

module load common
module load bedtools
module load UCSC
module load samtools

cd $WORKDIR

for NAME in `cut -f 1 sampleList.txt`; do
  BAMFILE=${NAME}.bam
  BEDGRAPHFILE=${NAME}.bedGraph
  BIGWIGFILE=${NAME}.bigWig

  if [ ! -r ${BEDGRAPHFILE}.gz ]; then
    READCOUNT=`samtools view -c $BAMFILE | awk '{ print 1e6/$1 }'`
    bedtools genomecov -scale $READCOUNT -bg -split -trackline -trackopts "name=\"$NAME\"" -ibam $BAMFILE -g $GENOMEFILE |grep -v random |grep -v Un|grep -v alt > $BEDGRAPHFILE
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

