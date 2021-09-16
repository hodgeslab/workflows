#!/usr/bin/env Rscript

# File to make peak heatmaps based on a matrix of read densities across peaks
# output from bwtool. If you use this code, please cite the following:
# Hodges et al., Nat Struct Mol Biol 25(1): 61-72 (2018), PMID: 29323272.
#
# To use, edit the following parameters below: conditions, tilesize, maxdist, normConst
# Whenever relevant, these values should match those used in the previous steps by bwtool.

library(hodgeslabR)
library(reshape)
library(ggplot2)

args <- commandArgs(TRUE)
cName <- args[1]
tName <- args[2]
sizefn <- args[3]
bedName <- args[4]

colorLen <- 64

options(bitmapType = "cairo")

# conditions <- c("treat-increased","treat-decreased") #, "unchanged")
tilesize <- 10 # for bwtool tile-averages
maxdist <- 4000

# define colormaps

# colormapControl <- colorRampPalette(c("#8bcde6","#c5dee6","#ffffff","#f6c54a","#a42018"))(colorLen)
# colormapTreat <- colorRampPalette(c("#8bcde6","#c5dee6","#ffffff","#f6c54a","#a42018"))(colorLen)

# colorControl <- "#a42018"
# colorTreat <- "#8bcde6"

colorControl <- "#808080"
colorTreat <- "#000000"

# read in size factors previously output from DESeq2 script
sizeFactors <- read.delim(sizefn,row.names=1,header=T)
sc1 <- sizeFactors$value[row.names(sizeFactors) == paste0(cName,"_rep1")]
sc2 <- sizeFactors$value[row.names(sizeFactors) == paste0(cName,"_rep2")]
st1 <- sizeFactors$value[row.names(sizeFactors) == paste0(tName,"_rep1")]
st2 <- sizeFactors$value[row.names(sizeFactors) == paste0(tName,"_rep2")]

# useful functions needed later
rotate <- function(x) t(apply(x, 2, rev))
procFile <- function(fn) {
  mat <- read.table(fn, header=F, row.names=NULL, sep="\t", fill=T)
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat[ ,-(1:6)])
}

# control1
fn <- paste0(bedName,"_",cName,"_rep1.txt")
imgControl1 <- procFile(fn)

# control2
fn <- paste0(bedName,"_",cName,"_rep2.txt")
imgControl2 <- procFile(fn)

# treat1
fn <- paste0(bedName,"_",tName,"_rep1.txt")
imgTreat1 <- procFile(fn)

# treat2
fn <- paste0(bedName,"_",tName,"_rep2.txt")
imgTreat2 <- procFile(fn)

imgControl <- (imgControl1/sc1 + imgControl2/sc2)/2
imgTreat <- (imgTreat1/st1 + imgTreat2/st2)/2

profControl <- colMeans(imgControl)
profTreat <- colMeans(imgTreat)
pos <- tilesize * (1:ncol(imgControl) - round(ncol(imgControl)/2))

# normalize by the mean value maxdist from center (values are "over background")
profControl <- profControl * 2 / (profControl[1] + profControl[length(profControl)])
profTreat <- profTreat * 2 / (profTreat[1] + profTreat[length(profTreat)])

dfControl <- data.frame(pos=pos,value=profControl,condition="control")
dfTreat <- data.frame(pos=pos,value=profTreat,condition="treat")
df <- rbind(dfControl, dfTreat)

p1 <- ggplot(df,aes(x=pos/1000,y=value,color=condition)) +
  theme_hodgeslab_basic() +
  geom_line(size=0.5*linescale) +
  scale_color_manual(values=c(colorControl, colorTreat), guide=F) +
  xlab("Distance from\nmotif (kb)") +
  ylab("Mean signal\n(over background)") +
  scale_x_continuous(breaks=c(-maxdist,0,maxdist)/1000) +
  facet_wrap(~ condition, nrow=1) +
  theme(panel.spacing = unit(0.5, "lines"))

pdf(paste0(bedName,"_metaGene.pdf"), height=0.55*2, width=0.6*2, useDingbats=F)
print(p1)
dev.off()

