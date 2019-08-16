#!/usr/bin/env Rscript

library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.NCBI.GRCh38)

args <- commandArgs(TRUE)
cName <- args[1]
tName <- args[2]

controlBams <- list.files(path="..", pattern=paste0(cName,"_rep[0-9].bam"))
treatBams <- list.files(path="..", pattern=paste0(tName,"_rep[0-9].bam"))
bamFiles <- paste0("../",c(controlBams,treatBams))

# load peaks, and use hack to remove "chr" prefix to make compatible with NCBI style chr names
peakfile <- "mergedPeaks.bed"
peaks <- getPeaks(peakfile, sort_peaks=T)
seqlevels(peaks) <- gsub("chr","",seqlevels(peaks))
peaks <- resize(peaks, width=300, fix="center")

fragment_counts <- getCounts(bamFiles, peaks,
                             paired =  TRUE,
                             by_rg = F,
                             format = "bam",
                             colData = DataFrame(celltype = c(rep("control",length(controlBams)),rep("treat",length(treatBams)))))

fragment_counts <- addGCBias(fragment_counts, 
                             genome = BSgenome.Hsapiens.NCBI.GRCh38)
# counts_filtered <- filterSamples(fragment_counts, min_depth = 1500,
#                                  min_in_peaks = 0.15)
counts_filtered <- filterPeaks(fragment_counts)
motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, counts_filtered,
                        genome = BSgenome.Hsapiens.NCBI.GRCh38)

# computing deviations
dev <- computeDeviations(object = counts_filtered, 
                         annotations = motif_ix)

res <- deviations(dev)
write.table(res,"chromVAR_deviations.txt",quote=F,sep="\t")

res <- deviationScores(dev)
write.table(res,"chromVAR_deviationScores.txt",quote=F,sep="\t")

