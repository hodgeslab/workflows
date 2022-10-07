#!/bin/sh

# For each R script below, the final three color definitions are optional.
# Optional colors for each heatmap correspond to low, medium, and high intensity.
# The BEDFILE variable corresponds to the bed filename used in sampleList.txt.

# Differential plots expect two replicates with suffixes _rep1 and _rep2.
# These plots also require a tab-delimited file containing the size factors
# from DESeq2 with the format below:
# samplename	sizefactor

# Annotation plots do not require a size factor and the scripts do not expect any
# suffix; be sure to use the full sample name (e.g. include _rep1 if needed).

module load R

# The entire content below this line (including R script calls) can be repeated
# for each bed file.
BEDFILE=sites.bed
BEDNAME=$(basename $BEDFILE .bed)

# Differential plots (conditions are quantitatively coupled for comparison):
Rscript plot_heatmaps_differential.R ATAC_THP1_DMSO ATAC_THP1_BRM014 sizeFactors.txt $BEDNAME
Rscript plot_heatmaps_differential.R ChIP_PU1_THP1_DMSO ChIP_PU1_THP1_BRM014 sizeFactors.txt $BEDNAME "#ffffff" "#e0e0e0" "#f50f02"

# Single plots (for annotation only):
Rscript plot_heatmap_single.R ChIP_SMARCA4_MOLM13_rep1 $BEDNAME "#ffffff" "#e0e0e0" "#0341fc"
Rscript plot_heatmap_single.R ChIP_H3K4me1_THP1_rep1 $BEDNAME "#ffffff" "#e0e0e0" "#02f740"
Rscript plot_heatmap_single.R ChIP_H3K4me3_THP1_rep1 $BEDNAME "#ffffff" "#e0e0e0" "#02f740"
Rscript plot_heatmap_single.R ChIP_H3K27ac_THP1_rep1 $BEDNAME "#ffffff" "#e0e0e0" "#02f740"
