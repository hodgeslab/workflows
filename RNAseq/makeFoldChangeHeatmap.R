#!~/usr/bin/Rscript

dataFile <- "diffexpr-results.txt"

library(gplots)
#library(RColorBrewer)
#my_palette <- colorRampPalette(brewer.pal(255,'RdBu'))(n = 255)
#my_palette <- rev(my_palette)
my_palette <- colorRampPalette(c("#0071bc", "gray95", "#d85218"))(n = 255)

# load data and set factor order based on order of appearance
df1 <- read.table(dataFile, header=T, row.names=NULL, sep="\t")
df1 <- na.omit(df1)
# df1 <- df1[df1$padj < 0.05, ]

# number of datasets
numDatasets <- 4

dataTable <- as.matrix(df1[ ,(dim(df1)[2]-numDatasets+1):dim(df1)[2]])
dataTable <- log2(dataTable / rowMeans(dataTable))
dataTable[!is.finite(dataTable)] <- NA
dataTable <- na.omit(dataTable)


clust <- kmeans(dataTable,3)
dataTable <- dataTable[order(clust$cluster), ]
labels <- df1$label[order(clust$cluster)]

colnames(dataTable) <- gsub("_reads","",colnames(dataTable))

# define scaling factor to get appropriate line widths in PDFs
# at present, defined as converting from mm into "R line widths" (1/96 of an inch):
# (25.4 mm/72.27 bigpts)*(96 linepts/72 pts) = 96/72/ggplot2:::.pt
linescale <- 96/72/ggplot2:::.pt

source('~/Projects/esBAF/Rlib/ggplot2_themes/themes_hodges.R')

pdf("outputTemplate.pdf",useDingbats=F,height=5,width=4)
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=0.75,cex.axis=0.75)


heatmap.2(dataTable, col=my_palette, margins=c(7,5), labRow=NA, cexCol=0.75, srtCol=45,
          density.info="none", trace="none", key=T, key.title="NA", key.xlab="log2 over mean",
          dendrogram="none", symkey=T, scale="none", Rowv=F, Colv=F, breaks=seq(-1,1,length.out=256))
#key.xtickfun=function() { return(list(at=c(-0,0.5,1),labels=as.character(c(-1,0,1)))) })
title(main=list(paste0(dim(dataTable)[1]," selected genes"), cex=0.75))

dev.off()

#ggsave("outputTemplate.pdf",useDingbats=F,height=1.5,width=2)



