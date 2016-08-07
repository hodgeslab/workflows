#!/bin/sh

EXPERIMENTS="All_conditions"

for i in $EXPERIMENTS; do
# rm -rf $i/merged*
# sh 3_makeDESeq2Table.sh $i.txt $i
 cd $i; Rscript ~/templates/DESeq2/ChIP/compare_control_treat_heatscatter.R rawCounts.txt 2; cd ..;
done


# make peak calls based on DESeq2 output

PTHRESH=0.10
FTHRESH=0.58496250072

for DATASET in */diffexpr-results.txt; do
  DATASETNAME=`echo $DATASET | cut -f 1 -d '/'`
  echo $DATASETNAME

  cat $DATASET | awk -v pthresh=$PTHRESH -v fthresh=$FTHRESH 'BEGIN{OFS="\t"; FS="\t" } NR > 1 { if($7 < pthresh && $3 <= -fthresh) print $9,$10,$11,$8,$7,$3 }' > $DATASETNAME/treat-decreased.bed
  cat $DATASET | awk -v pthresh=$PTHRESH -v fthresh=$FTHRESH 'BEGIN{OFS="\t"; FS="\t" } NR > 1 { if($7 < pthresh && $3 >= fthresh) print $9,$10,$11,$8,$7,$3 }' > $DATASETNAME/treat-increased.bed
  cat $DATASET | awk -v pthresh=$PTHRESH -v fthresh=$FTHRESH 'BEGIN{OFS="\t"; FS="\t" } NR > 1 { if($7 > pthresh || $7 == "NA" || ($3 > -fthresh && $3 < fthresh)) print $9,$10,$11,$8,$7,$3 }' > $DATASETNAME/unchanged.bed

done
