#!/usr/bin/env Rscript

# File to make peak heatmaps based on a matrix of read densities across peaks
# output from bwtool. If you use this code, please cite the following:
# Hodges et al., Nat Struct Mol Biol 25(1): 61-72 (2018), PMID: 29323272.
#
# To use, edit the following parameters below: conditions, tilesize, maxdist, normConst
# Whenever relevant, these values should match those used in the previous steps by bwtool.

library(RColorBrewer)
library(reshape)

args <- commandArgs(TRUE)
cName <- args[1]
tName <- args[2]
sizefn <- args[3]

options(bitmapType = "cairo")

conditions <- c("treat-increased","treat-decreased") #, "unchanged")
tilesize <- 10 # for bwtool tile-averages
maxdist <- 4000
colorLen <- 64

# define colormaps
# colormapControl <- brewer.pal("Reds",n=9)
# colormapControl[1] <- "#FFFFFF"
# colormapControl[10] <- colormapControl[9]
# colormapTreat <- c("#FFFFFF",brewer.pal("Blues",n=9))
# colormapTreat[1] <- "#FFFFFF"
# colormapTreat[10] <- colormapTreat[9]

c1 <- "#ffffff"
c2 <- "#e0e0e0"
c3 <- "#000000"

colormapControl <- colorRampPalette(c(c1,c2,c3))(colorLen)
colormapTreat <- colorRampPalette(c(c1,c2,c3))(colorLen)

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

for(i in conditions) {
  
  # control1
  fn <- paste0(i,"_",cName,"_rep1.txt")
  imgControl1 <- procFile(fn)

  # control2
  fn <- paste0(i,"_",cName,"_rep2.txt")
  imgControl2 <- procFile(fn)

  # treat1
  fn <- paste0(i,"_",tName,"_rep1.txt")
  imgTreat1 <- procFile(fn)

  # treat2
  fn <- paste0(i,"_",tName,"_rep2.txt")
  imgTreat2 <- procFile(fn)

  imgControl <- (imgControl1/sc1 + imgControl2/sc2)/2
  imgTreat <- (imgTreat1/st1 + imgTreat2/st2)/2
  
  # identify peak: within middle 50%, sum both control and treat, find peak, and recenter
  startPos <- t(apply(imgControl+imgTreat,1,function(x) {
    openRange = c(round(0.25*length(x)),round(0.75*length(x)))
    peaki = which.max(x[openRange[1]:openRange[2]])
  }))
  
  imgControl2 <- matrix(NA,nrow = nrow(imgControl),ncol = round(0.5*ncol(imgControl)))
  imgTreat2 <- matrix(NA,nrow = nrow(imgTreat),ncol = round(0.5*ncol(imgTreat)))
  for(j in 1:nrow(imgControl)) {
    imgControl2[j,] <- imgControl[j, startPos[j]:(startPos[j]+round(0.5*ncol(imgControl))-1)]
    imgTreat2[j,] <- imgTreat[j, startPos[j]:(startPos[j]+round(0.5*ncol(imgControl))-1)]
  }
  
  imgControl <- imgControl2
  imgTreat <- imgTreat2
  
  # Set a normalization constant. This is a way of controlling image brightness
  # and contrast (Z scale).
  imgTot <- cbind(imgControl,imgTreat)
  imgTot <- t(scale(t(imgTot)))
  imgControl <- imgTot[,1:ncol(imgControl)]
  imgTreat <- imgTot[,(ncol(imgControl)+1):(ncol(imgControl)+ncol(imgTreat))]
  
  # sort (based on middle intensity)
  sel <- c(round(0.45*ncol(imgControl)),round(0.55*ncol(imgControl)))
  idx <- order(pmax(rowSums(imgControl[,sel],rowSums(imgTreat[ ,sel]))),decreasing = T)
  # imgControl <- imgControl[idx, ]
  # imgTreat <- imgTreat[idx, ]
  
  # cap maximum intensities
  normConst <- 4
  imgControl[imgControl > normConst] <- normConst
  imgTreat[imgTreat > normConst] <- normConst
  
  # write images
  outfile <- paste0(i,'_control.tif')                # 100 pixels wide, and
  tiff(outfile,w=100,h=round(dim(imgControl)[1]/1)) # each verticle pixel represents 1 sites
  par(mar = c(0,0,0,0))
  image(rotate(as.matrix(imgControl)), col=colormapControl,axes = FALSE,breaks=c(seq(0,normConst,length.out=colorLen),100))
  dev.off()
  
  outfile <- paste0(i,'_treat.tif')                # 100 pixels wide, and
  tiff(outfile,w=100,h=round(dim(imgTreat)[1]/1)) # each verticle pixel represents 1 sites
  par(mar = c(0,0,0,0))
  image(rotate(as.matrix(imgTreat)), col=colormapTreat,axes = FALSE,breaks=c(seq(0,normConst,length.out=colorLen),100))
  dev.off()
}
