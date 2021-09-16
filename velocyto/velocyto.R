#!~/usr/bin/Rscript

library(methods)
library(hodgeslabR)
library(velocyto.R)
library(pagoda2)
library(igraph)
library(pheatmap)
library(grid)
library(RColorBrewer)

args <- commandArgs(TRUE)
infile <- args[1]

ldat <- read.loom.matrices(infile)
emat <- ldat$spliced
# filter on cells with large numbers of splicing events
emat <- emat[,colSums(emat)>=1e3]

# eliminate duplicated gene names (this is mostly for snoRNA and Y_RNA genes)
emat <- emat[!duplicated(rownames(emat)), ]

r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)

pdf('1_overdisperse_genes.pdf', useDingbats=F, height=3, width=6)
r$adjustVariance(plot=T,do.par=T,gam.k=10)
dev.off()

r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)

r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)

pdf('2_clusters_and_proliferating.pdf', useDingbats=F, height=3, width=6)
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"AURKA"],main='AURKA')
dev.off()

# VELOCITY ESTIMATION

emat <- ldat$spliced; nmat <- ldat$unspliced
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]; # restrict to cells that passed p2 filter
# take cluster labels
cluster.label <- r$clusters$PCA[[1]]
cell.colors <- pagoda2:::fac2col(cluster.label)
# take embedding
emb <- r$embeddings$PCA$tSNE

# In addition to clustering and the t-SNE embedding, from the p2 processing we will also take a cell-cell
# distance, which will be better than the default whole-transcriptome correlation distance that velocyto.R
# would normally use. The 1-correlation here just turns correlation matrix into a dissimilarity score:
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

lims <- c(-0.6,0.6)
mat <- as.matrix(1-cell.dist)
diag(mat) <- 1 # diagonal is perfectly correlated
mat[mat < lims[1]] <- lims[1]
mat[mat > lims[2]] <- lims[2]
p1 <- pheatmap(mat, show_rownames=F, show_colnames=F, silent=T,
  color = rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100)),
  breaks = seq(lims[1],lims[2],length.out = 101), legend_breaks = c(-0.5,0,0.5))
#  annotation_row = data.frame(cluster = cluster.label),
#  annotation_col = data.frame(cluster = cluster.label))
p1$gtable$grobs[[1]]$gp <- gpar(lwd = 0.5*linescale)
p1$gtable$grobs[[2]]$gp <- gpar(lwd = 0.5*linescale)

# pdf('3_cell_dist.pdf', useDingbats=F, height=3, width=3.5)
# grid.newpage()
# grid.draw(p1$gtable)
# dev.off()

png('3_cell_corr.png', height=3, width=3.5, units='in', res=300)
grid.newpage()
grid.draw(p1$gtable)
dev.off()

# Filter genes based on the minimum average expresion magnitude (in at least one of the clusters),
# output total number of resulting valid genes:
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

# Estimate RNA velocity (using gene-relative model with k=20 cell kNN pooling and using
# top/bottom 2% quantiles for gamma fit):
fit.quantile <- 0.02
# rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)
# k=100 cell neighborhood
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=100,cell.dist=cell.dist,fit.quantile=fit.quantile)

pdf('4_velocity_field.pdf', useDingbats=F, height=5, width=5)
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),
  cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()

# Visualize a fit for a particular gene (we reuse rvel.cd to save on calcualtions here):
gene <- "FN1"
pdf('5_gene_fit.pdf', useDingbats=F, height=2, width=7)
gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 100,kGenes=1,fit.quantile=fit.quantile,cell.emb=emb,
  cell.colors=cell.colors,cell.dist=cell.dist,show.gene=gene,old.fit=rvel.cd,do.par=T)
dev.off()
