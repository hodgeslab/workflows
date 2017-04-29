#!/bin/sh

WORKDIR=XXWORKDIRXX
GENOME=mm9

cd $WORKDIR
mkdir -p scripts

for FN in *.fastq.gz; do
  NAME=`basename $FN | sed s/.fastq.gz//`

  cat > scripts/map_${NAME}.qsub <<EOL
WORKDIR=$WORKDIR
cd \$WORKDIR
if [ -r ${FN} ]; then gunzip ${FN}; fi
if [ ! -r ${NAME}.map ]; then ~/bin/bowtie -m 1 ${GENOME} ${NAME}.fastq ${NAME}.map; fi
cat ${NAME}.map | awk 'BEGIN{FS=OFS="\t"}{if(NR==1) for(i=1; i<=NF; i++) if(\$i~/^chr/) chr=i; start=\$(chr+1); len=length(\$(chr+2)); end=(\$(chr+1)+len); if(start<end) print \$chr,start,end,len,\$1,\$(chr-1) }' | \\
  sort -k1,1 -k2,2n -k3,3n -k6,6 -u > ${NAME}_sort_noDups.bed
EOL
done
