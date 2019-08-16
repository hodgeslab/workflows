#!/bin/sh

INFILE="All_conditions.txt"
GENOME=hg38
WORKDIR=`pwd`
TEMPLATE=slurm_DESeq2

module load R

function all_vs_all () {
  iseen=""
  for irep in `cat $INFILE`; do
    iname=`sed s/_rep[0-9]*// <<< $irep`

    if grep -q "$iname" <<< $iseen; then
      continue;
    fi

    jseen=""
    for jrep in `cat $INFILE`; do
      jname=`sed s/_rep[0-9]*// <<< $jrep`
      if ! (grep -q "$jname" <<< $iseen || grep -q "$jname" <<< $iname || grep -q "$jname" <<< $jseen); then
        exptname="${iname}-${jname}";
        # assemble txt file from i and j
        echo "$exptname"
        grep $iname $INFILE > ${exptname}.txt
        grep $jname $INFILE >> ${exptname}.txt
      fi
    jseen="$jseen $jname"
    done
  iseen="$iseen $iname";
  done
}

# make all vs. all comparisons
EXPERIMENTS="All_conditions `all_vs_all`"

for EXPERIMENT in $EXPERIMENTS; do
  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXEXPERIMENTXX%$EXPERIMENT% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/DESeq2_${EXPERIMENT}.slurm
done
