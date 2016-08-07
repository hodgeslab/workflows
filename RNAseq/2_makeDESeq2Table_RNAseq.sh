#!/bin/sh

PATH=$PATH:~/src/bedtools2/bin
WORKDIR=.
INPUTDIR=.
SAMPLELIST=$1
EXPERIMENTNAME=$2
OUTPUTDIR=$EXPERIMENTNAME
GENOME=mm9

cd $WORKDIR

mkdir -p $OUTPUTDIR

for SAMPLE in `cat $SAMPLELIST`; do
  echo "Processing $SAMPLE ..."
  cp $INPUTDIR/${SAMPLE}_counts.txt $OUTPUTDIR
done

# partially assemble tableFile

echo "Assembling data from all samples ..."

TABLEFILE=$OUTPUTDIR/sampleTable.txt
GENEBODIESFILE=$OUTPUTDIR/geneBodies.txt

echo -e "sampleName\tfileName" > $TABLEFILE
for SAMPLE in `cat $SAMPLELIST`; do
  echo -e "${SAMPLE}\t${SAMPLE}_counts.txt" >> $TABLEFILE
done

cat ~/data/UCSC_Downloads/${GENOME}/refGene_${GENOME}.gtf | grep "refGene	transcript"|sed 's/gene_id\ \"\([^\"]*\)\".*/\1/'  | sort -k 9,9 | \
awk 'BEGIN{FS=OFS="\t"; i=0}{if(!seen[$9]) {i++; gene[i] = $9; chr[i]=$1; start[i]=$4; end[i]=$5; strand[i]=$7; seen[$9]++} else {if($4 < start[i]) start[i]=$4; if($5 > end[i]) end[i]=$5}} END{print "id\tchr\tstart\tend\tstrand"; for(j=1;j<=i;j++) print gene[j],chr[j],start[j],end[j],strand[j] }' > $GENEBODIESFILE
