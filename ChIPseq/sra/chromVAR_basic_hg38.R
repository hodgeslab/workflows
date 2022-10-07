#!/usr/bin/env Rscript

library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(hodgeslabR)
library(ggplot2)
library(ggrepel)

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

df <- as.data.frame(res)
csel <- 1:length(controlBams)
tsel <- (1+length(controlBams)):ncol(df)
df$controlmean <- rowMeans(df[,csel])
df$treatmean <- rowMeans(df[,tsel])
df$motif <- row.names(df)
df$deltaScoreMean <- df$treatmean - df$controlmean
df$deltaScore1 <- df[,tsel[1]] - df[,csel[1]]
df$deltaScore2 <- df[,tsel[2]] - df[,csel[2]]
df$rank <- rank(df$deltaScoreMean,ties.method = "first")
df <- df[order(df$rank,decreasing = F), ]

ymax <- max(abs(df$deltaScoreMean))
df <- melt(df,id.vars = c("rank","motif"), measure.vars = c("deltaScore1","deltaScore2","deltaScoreMean"))

df$color <- "#000000"
# assign other colors, e.g.
# sel <- grepl("FOS",df$motif))
# df$color[sel] <- "#fd902a"

df$alpha <- 1
df$alpha[df$variable == "deltaScore1"] <- 0.2
df$alpha[df$variable == "deltaScore2"] <- 0.2

df$shape <- 16
df$shape[df$variable == "deltaScoreMean"] <- 18

# reorder plot (front color to back black)
df <- rbind(df[df$color=="#000000",],df[df$color != "#000000",])
df$size <- 1
df$size[df$color != "#000000"] <- 3

df$label <- df$motif
df$label[df$rank > 10 & df$rank < max(df$rank) - 10 ] <- ""
df$label[df$variable != "deltaScoreMean"] <- ""

p1 <- ggplot(df,aes(x=rank,y=value,color=color,alpha=alpha,label=label)) +
  theme_hodgeslab_basic() +
  geom_hline(aes(yintercept = 0),color="#000000",linetype="dotted") +
  geom_point(aes(size=df$size,stroke=0,shape=df$shape)) +
  scale_size_identity() +
  scale_shape_identity() +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_text_repel(size = 1, segment.size = 0.25*linescale,
                  direction="y") +
  xlab("Motif rank") +
  ylab("ATAC score diff") +
  ylim(-1.2*ymax,1.2*ymax)

pdf("chromVAR_deviationScores.pdf",useDingbats = F,height = 1.5,width=1.75)
print(p1)
dev.off()
