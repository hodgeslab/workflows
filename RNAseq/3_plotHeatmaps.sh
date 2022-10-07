#!/bin/sh

EXPERIMENTS="All_conditions"
WORKDIR=$(pwd)
SCRIPTFILE=plotDiffHeatmap.R

module load R

cd $WORKDIR
for i in $EXPERIMENTS; do

 CONTROLSTR=$(head -1 ${i}.txt | sed s/'_rep1'//)
 CONTROLNUM=$(grep "$CONTROLSTR" ${i}.txt | wc -l)

 cd $i; Rscript ../${SCRIPTFILE}; cd ..;
done

tar cvf data.tar ${EXPERIMENTS}

