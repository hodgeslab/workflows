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
c_name <- args[1]
t_name <- args[2]
size_fn <- args[3]
bed_name <- args[4]

# if colors are not defined, use defaults
if(length(args) == 7) {
  c1 <- args[5]
  c2 <- args[6]
  c3 <- args[7]
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

# read in size factors previously output from DESeq2 script
size_factors <- read.delim(size_fn,row.names=1,header=FALSE)
colnames(size_factors) <- "value"

size_c1 <- size_factors$value[row.names(size_factors) == paste0(c_name,"_rep1")]
size_c2 <- size_factors$value[row.names(size_factors) == paste0(c_name,"_rep2")]
size_t1 <- size_factors$value[row.names(size_factors) == paste0(t_name,"_rep1")]
size_t2 <- size_factors$value[row.names(size_factors) == paste0(t_name,"_rep2")]

# useful functions needed later
rotate <- function(x) t(apply(x, 2, rev))
proc_file <- function(fn) {
  mat <- read.table(fn, header=F, row.names=NULL, sep="\t")
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat[ ,-(1:6)])
}

# control 1
fn <- paste0(bed_name,"_",c_name,"_rep1.txt")
img_c1 <- proc_file(fn)

# control 2
fn <- paste0(bed_name,"_",c_name,"_rep2.txt")
img_c2 <- proc_file(fn)

# treat 1
fn <- paste0(bed_name,"_",t_name,"_rep1.txt")
img_t1 <- proc_file(fn)

# treat 2
fn <- paste0(bed_name,"_",t_name,"_rep2.txt")
img_t2 <- proc_file(fn)

# average across both replicates, adjusted by size factor
img_c <- (img_c1/size_c1 + img_c2/size_c2)/2
img_t <- (img_t1/size_t1 + img_t2/size_t2)/2
  
# identify peak: within middle 50%, sum both control and treat, find peak, and recenter
start_pos <- t(apply(img_c + img_t,1,function(x) {
  open_range = c(round(0.25*length(x)),round(0.75*length(x)))
  peaki = which.max(x[open_range[1]:open_range[2]])
}))

img_c2 <- matrix(NA,nrow = nrow(img_c),ncol = round(0.5*ncol(img_c)))
img_t2 <- matrix(NA,nrow = nrow(img_t),ncol = round(0.5*ncol(img_t)))

for(j in 1:nrow(img_c)) {
  img_c2[j,] <- img_c[j, start_pos[j]:(start_pos[j]+round(0.5*ncol(img_c))-1)]
  img_t2[j,] <- img_t[j, start_pos[j]:(start_pos[j]+round(0.5*ncol(img_t))-1)]
}
  
img_c <- img_c2
img_t <- img_t2

img_tot <- cbind(img_c, img_t)
low_lim <- apply(img_tot,1,function(x) { as.numeric(quantile(x,0.15)) } )
high_lim <- apply(img_tot,1,function(x) { as.numeric(quantile(x,0.99)) } )
img_tot <- (img_tot-low_lim)/(high_lim-low_lim)
img_c <- img_tot[,1:ncol(img_c)]
img_t <- img_tot[,(ncol(img_c)+1):(ncol(img_c)+ncol(img_t))]

# write images
out_file <- paste0(bed_name,'_', c_name, '.tif') # 100 pixels wide (spanning max_dist), and
tiff(out_file,w=100,h=round(dim(img_c)[1]/1))    # each verticle pixel represents 1 site
par(mar = c(0,0,0,0))
image(rotate(as.matrix(img_c)), col=colormap,axes = FALSE,breaks=c(seq(0,1,length.out=color_len),100))
dev.off()

out_file <- paste0(bed_name,'_', t_name, '.tif') # 100 pixels wide (spanning max_dist), and
tiff(out_file,w=100,h=round(dim(img_t)[1]/1))    # each verticle pixel represents 1 site
par(mar = c(0,0,0,0))
image(rotate(as.matrix(img_t)), col=colormap,axes = FALSE,breaks=c(seq(0,1,length.out=color_len),100))
dev.off()

