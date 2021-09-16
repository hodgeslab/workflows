#!/bin/sh

BIGWIG1=$1
BIGWIG2=$2
SAMPLENAME=$3
BEDFILE=$4

MAXDIST=4000
TILE=10

module load bedtools
module load common

makeMatrix() {
  BEDFILE=$1
  BIGWIG1=$2
  BIGWIG2=$3
  MAXDIST=$4
  TILE=$5
  SAMPLENAME=$6

  BEDNAME=$(basename $BEDFILE .bed)
  OUT1=${BEDNAME}_${SAMPLENAME}_rep1.txt
  OUT2=${BEDNAME}_${SAMPLENAME}_rep2.txt
  LOCALBED=local_`basename $BEDFILE`

  cut -f 1-3 $BEDFILE > $LOCALBED
  bwtool -keep-bed matrix ${MAXDIST}:${MAXDIST} $LOCALBED -tiled-averages=${TILE} $BIGWIG2 $OUT1
  sed -i 's/NA/0/g' $OUT1
  tr "," "\t" < $OUT1 > $OUT1.tmp
  mv $OUT1.tmp $OUT1

  bwtool -keep-bed matrix ${MAXDIST}:${MAXDIST} $LOCALBED -tiled-averages=${TILE} $BIGWIG1 $OUT2
  sed -i 's/NA/0/g' $OUT2
  tr "," "\t" < $OUT2 > $OUT2.tmp
  mv $OUT2.tmp $OUT2
}

makeMatrix $BEDFILE $BIGWIG1 $BIGWIG2 $MAXDIST $TILE $SAMPLENAME
