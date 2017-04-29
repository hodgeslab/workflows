#!/bin/sh

EXPERIMENT=$1
CONTROLBIGWIG=$2
TREATBIGWIG=$3
MAXDIST=4000
TILE=100

PATH=$PATH:~/src/bedtools2/bin

for NAME in treat-increased treat-decreased unchanged; do
  BEDFILE=${EXPERIMENT}/${NAME}.bed
  CONTROLOUTFILE=${EXPERIMENT}/${NAME}_control.txt
  TREATOUTFILE=${EXPERIMENT}/${NAME}_treat.txt
  CLEANBEDFILE=$(mktemp ${BEDFILE}.tmp.XXXXXX)

  cut -f 1-3 $BEDFILE > $CLEANBEDFILE

  bwtool matrix ${MAXDIST}:${MAXDIST} -tiled-averages=$TILE $CLEANBEDFILE $CONTROLBIGWIG $CONTROLOUTFILE
  bwtool matrix ${MAXDIST}:${MAXDIST} -tiled-averages=$TILE $CLEANBEDFILE $TREATBIGWIG $TREATOUTFILE
  rm "$CLEANBEDFILE"

done
