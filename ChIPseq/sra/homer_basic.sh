#!/bin/sh

GENOME=$1

module load homer

# HOMER
CONTROLSITES=unchanged.bed
TREATSITES=treat-increased.bed
OUTPUTDIR=homerOutput_treat-increased
findMotifsGenome.pl $TREATSITES $GENOME $OUTPUTDIR -bg $CONTROLSITES -size given
annotatePeaks.pl $TREATSITES $GENOME -go ${OUTPUTDIR} > $OUTPUTDIR/goStdout_treat.txt
GenomeOntologyControl.pl $TREATSITES $CONTROLSITES $GENOME $OUTPUTDIR

mkdir -p ${OUTPUTDIR}/control
annotatePeaks.pl $CONTROLSITES $GENOME -go ${OUTPUTDIR}/control > $OUTPUTDIR/control/goStdout_treat.txt

CONTROLSITES=unchanged.bed
TREATSITES=treat-decreased.bed
OUTPUTDIR=homerOutput_treat-decreased
findMotifsGenome.pl $TREATSITES $GENOME $OUTPUTDIR -bg $CONTROLSITES -size given
annotatePeaks.pl $TREATSITES $GENOME -go ${OUTPUTDIR} > $OUTPUTDIR/goStdout_treat.txt
GenomeOntologyControl.pl $TREATSITES $CONTROLSITES $GENOME $OUTPUTDIR

mkdir -p ${OUTPUTDIR}/control
annotatePeaks.pl $CONTROLSITES $GENOME -go ${OUTPUTDIR}/control > $OUTPUTDIR/control/goStdout_treat.txt

CONTROLSITES=treat-decreased.bed
TREATSITES=treat-increased.bed
OUTPUTDIR=homerOutput_treat-increased-decreased
findMotifsGenome.pl $TREATSITES $GENOME $OUTPUTDIR -bg $CONTROLSITES -size given
GenomeOntologyControl.pl $TREATSITES $CONTROLSITES $GENOME $OUTPUTDIR
