#!/usr/bin/env Rscript

# File to make peak heatmaps based on a matrix of read densities across peaks
# output from bwtool. If you use this code, please cite the following:
# Hodges et al., Nat Struct Mol Biol 25(1): 61-72 (2018), PMID: 29323272.
#
# Always ensure that tile_size and max_dist match those used in the matrix
# calculation performed by bwtool.

library(RColorBrewer)
library(reshape)

args <- commandArgs(TRUE)
sample_name <- args[1]
bed_name <- args[2]

# if colors are not defined, use defaults
if(length(args) == 5) {
  c1 <- args[3]
  c2 <- args[4]
  c3 <- args[5]
} else {
  c1 <- "#ffffff"
  c2 <- "#e0e0e0"
  c3 <- "#000000"
}

options(bitmapType = "cairo")

tile_size <- 10 # for bwtool tile-averages
max_dist <- 4000

# define colormap
color_len <- 64
colormap <- colorRampPalette(c(c1,c2,c3))(color_len)

# useful functions needed later
rotate <- function(x) t(apply(x, 2, rev))
proc_file <- function(fn) {
  mat <- read.table(fn, header=F, row.names=NULL, sep="\t")
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat[ ,-(1:6)])
}

# single sample
fn <- paste0(bed_name,"_",sample_name,".txt")
img_c <- proc_file(fn)

# identify peak: within middle 50%, sum both control and treat, find peak, and recenter
start_pos <- t(apply(img_c,1,function(x) {
  open_range = c(round(0.25*length(x)),round(0.75*length(x)))
  peaki = which.max(x[open_range[1]:open_range[2]])
}))

img_c2 <- matrix(NA,nrow = nrow(img_c),ncol = round(0.5*ncol(img_c)))

for(j in 1:nrow(img_c)) {
  img_c2[j,] <- img_c[j, start_pos[j]:(start_pos[j]+round(0.5*ncol(img_c))-1)]
}
  
img_c <- img_c2

low_lim <- apply(img_c,1,function(x) { as.numeric(quantile(x,0.15)) } )
high_lim <- apply(img_c,1,function(x) { as.numeric(quantile(x,0.99)) } )
img_c <- (img_c-low_lim)/(high_lim-low_lim)

# write image
out_file <- paste0(bed_name,'_', sample_name, '.tif') # 100 pixels wide (spanning max_dist), and
tiff(out_file,w=100,h=round(dim(img_c)[1]/1))    # each verticle pixel represents 1 site
par(mar = c(0,0,0,0))
image(rotate(as.matrix(img_c)), col=colormap,axes = FALSE,breaks=c(seq(0,1,length.out=color_len),100))
dev.off()
