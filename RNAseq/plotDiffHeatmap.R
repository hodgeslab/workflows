#!/usr/bin/env Rscript

library(gplots)
library(RColorBrewer)

inFile <- "diffexpr-results_rld.txt"
df <- read.delim(inFile, header=T, sep="\t", stringsAsFactors = F)

# heatmap
# df <- resdata_rld
df <- df[df$padj < 0.25, ]
df <- na.omit(df)
dataTable <- as.matrix(df[ ,13:dim(df)[2]])
head(dataTable)
#dataTable <- log2(dataTable / rowMeans(dataTable))
dataTable <- t(scale(t(dataTable)))
dataTable[!is.finite(dataTable)] <- NA
dataTable <- na.omit(dataTable)

clust <- kmeans(dataTable,2)
dataTable <- dataTable[order(clust$cluster), ]
labels <- df$label[order(clust$cluster)]

colnames(dataTable) <- gsub("_reads","",colnames(dataTable))

pdf("diffexpr-log2heatmap.pdf",useDingbats=F,height=5,width=4)
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=0.75,cex.axis=0.75)
my_palette <- colorRampPalette(c("#0071bc", "gray95", "#d85218"))(n = 255)
heatmap.2(dataTable, col=my_palette, margins=c(7,5), labRow=labels, cexCol=0.75, srtCol=45,
          density.info="none", trace="none", key=T, key.title="NA", key.xlab="Z score",
          dendrogram="none", symkey=T, scale="none", Rowv=F, Colv=F, breaks=seq(-1,1,length.out=256))
#key.xtickfun=function() { return(list(at=c(-0,0.5,1),labels=as.character(c(-1,0,1)))) })
title(main=list(paste0(dim(dataTable)[1]," selected genes"), cex=0.75))

dev.off()
