#!~/usr/bin/Rscript

# assumes will read in tab-delimited file with one header line, with columns as the following:
# sampleName
# fileName

# command line args
args <- commandArgs(TRUE)
tableFile <- args[1]
numControl <- as.numeric(args[2])
# numTreat <- as.numeric(args[3])

# organize tabulated data
geneBodiesFile <- "geneBodies.txt"
geneBodies <- read.table(geneBodiesFile, header=T, row.names="id", sep="\t", stringsAsFactors=F)
sampleTable <- read.table(tableFile, header=T, row.names=NULL, sep="\t")
sampleTable$dataset <- sampleTable$sampleName

# populate conditions
numTreat <- nrow(sampleTable) - numControl
condition <- factor(c(rep("control", numControl), rep("treat", numTreat)))
sampleTable$condition <- condition

# load DESeq2 and process data
library(DESeq2)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".", design = ~condition) 
dds$condition <- relevel(dds$condition, "control")

# do everything automagically
# dds <- DESeq(dds, fitType="local")
# or manally
# because ChIP seq has background peak calls, we normalize not to housekeeping genes:
# dds <- estimateSizeFactors(dds)
# dds <- estimateSizeFactors(dds, controlGenes = (rowMeans(countdata) > 1000))
# but instead set the sizeFactors to be the total number of reads within peaks
# sizeFactors(dds) <- colSums(countdata) / 1e6
sizeFactors(dds) <- colSums(counts(dds)) / 1e6

# we can also ignore the bottom 10% of total peaks and assume them to be background
# bgFrac <- 0.10
# pseudoPeakOrder <- order(rowSums(colSums(countdata) * countdata / 1e6), decreasing = TRUE)
# pseudoPeakOrder <- pseudoPeakOrder[c(1:round((1 - bgFrac) * length(pseudoPeakOrder)))]
# sizeFactors(dds) <- colSums(countdata[pseudoPeakOrder, ]) / 1e6

# or we can only consider the top 5% of peaks to stay saturated
# topFrac <- 0.05
# pseudoPeakOrder <- order(rowSums(colSums(countdata) * countdata / 1e6), decreasing = TRUE)
# pseudoPeakOrder <- pseudoPeakOrder[c(1:round(topFrac * length(pseudoPeakOrder)))]
# sizeFactors(dds) <- colSums(countdata[pseudoPeakOrder, ]) / 1e6

dds <- estimateDispersions(dds, fitType = "local")
dds <- nbinomWaldTest(dds)

# plot dispersions
pdf("qc-dispersions.pdf")
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1.75,cex.axis=1.75)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

rld <- rlogTransformation(dds, blind=FALSE)
#rlddata <- as.data.frame(assay(rld), row.names=c(1:nrow(countdata)))
head(assay(rld))
# hist(assay(rld))
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# plot heatmap distance of samples
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
pdf("qc-heatmap-samples.pdf")
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1,cex.axis=1)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margins=c(12, 12), main="Sample Distance Matrix", cexRow=1,cexCol=1,srtCol=45)
dev.off()

# plot PCA
pdf("qc-pca.pdf", useDingbats=F, height=4,width=8)
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1.75,cex.axis=1.75)
DESeq2::plotPCA(rld, intgroup=c("condition", "dataset"))
dev.off()

res <- results(dds)
res$label <- row.names(res)
#res$chr <- dataChr
#res$startPos <- dataStart
#res$endPos <- dataEnd
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with gene body data
resdata <- merge(as.data.frame(res), geneBodies, by="row.names", sort=F)
resdata <- transform(resdata,row.names=Row.names, Row.names=NULL)
head(resdata)
## Merge with normalized count data
resdata <- merge(resdata, as.data.frame(counts(dds, normalized=T)), by="row.names", sort=F)
#resdata <- merge(as.data.frame(res), as.data.frame(rlddata), by="row.names", sort=F)
names(resdata)[1] <- "id"
head(resdata)

write.table(resdata, file="diffexpr-results.txt", sep="\t", quote=FALSE, row.names=F)

## output MA plot
#maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
#  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
#  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
#  if (labelsig) {
#    require(calibrate)
#    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=id, cex=textcx, col=2))
#  }
#}
pdf("diffexpr-maplot.pdf")
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1.75,cex.axis=1.75)
DESeq2::plotMA(res, ylim=c(-2,2), ylab="log2 fold change", xlab="mean expression (RPM)")
abline(h=c(-log2(1.5),log2(1.5)),col="dodgerblue",lwd=2)
#axis(side=1,at=c(10,100,1000))
dev.off()

# output volcanoplot
volcanoplot <- function (res, lfcthresh=1, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=2, ...) {
  par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1.75,cex.axis=1.75)
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, xlab="log2 fold change", ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=label, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
pdf("diffexpr-volcanoplot.pdf")
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8) # , xlim=c(-5, 5), ylim=c(-1, 55))
dev.off()

library(boot)

corbootf <- function(d,i) {
	d2 <- d[i, ]
	return(cor(d2$rep1, d2$rep2, method="pearson"))
}
plotFoldChangeCorrelation <- function (resdata, axisLim=1,sigthresh=0.10, numControl=2, numTreat=2, main=NULL, logthresh=log2(1.5)) {
	df <- as.data.frame(resdata[ ,c(7,13,14,13+numControl,14+numControl)])
	colnames(df) <- c("padj","control1","control2","treat1","treat2")
	df$rep1 <- log2(df$treat1+0.5) - log2(df$control1+0.5)
	df$rep2 <- log2(df$treat2+0.5) - log2(df$control2+0.5)
	plot(0,0,type="n",axes=FALSE,
		xlab="",xaxs="i",xlim=c(-axisLim,axisLim),
		ylab="",yaxs="i",ylim=c(-axisLim,axisLim))
	title(main=main,xlab="log2 fold change, replicate 1",ylab="log2 fold change, replicate 2")
	axis(1, at = pretty(c(-axisLim,axisLim)))
	axis(2, at = pretty(c(-axisLim,axisLim)))
	abline(h=0,v=0,col="black",lty=3)
	# with(subset(df, padj<sigthresh), points(rep1, rep2, pch=20, col=rgb(1,0,0,0.5)))
	# with(subset(df, padj>sigthresh), points(rep1, rep2, pch=20, col=rgb(0,0,0,0.075)))
	sigPoints <- df$padj < sigthresh & abs(resdata$log2FoldChange) > logthresh
	sigPointsIncr <- df$padj < sigthresh & resdata$log2FoldChange > logthresh
	sigPointsDecr <- df$padj < sigthresh & resdata$log2FoldChange < -logthresh
	nsigPoints <- ! sigPoints | is.na(df$padj)
#	points(df$rep1[sigPoints], df$rep2[sigPoints], pch=20, col=rgb(1,0,0,0.05))
#	points(df$rep1[nsigPoints], df$rep2[nsigPoints], pch=20, col=rgb(0,0,0,0.01))
	points(df$rep1[sigPointsIncr], df$rep2[sigPointsIncr], pch=20, col=rgb(0.8500,0.3250,0.0980,0.05))
	points(df$rep1[sigPointsDecr], df$rep2[sigPointsDecr], pch=20, col=rgb(0,0.4470,0.7410,0.05))
	points(df$rep1[nsigPoints], df$rep2[nsigPoints], pch=20, col=rgb(0,0,0,0.01))
	rval = cor(df$rep1, df$rep2, method="pearson")
	
	# bootstrapping with 100 replications 
	results <- boot(data=df, statistic=corbootf, R=100)
	rsqci <- boot.ci(results, type="basic")
	label <- sprintf("R = %.3f\n(%.3f, %.3f) 95%% CI",rval,rsqci$basic[4],rsqci$basic[5],3)
	text(-0.5*axisLim,0.75*axisLim,label,offset=0,cex=1.4)
	}

if(numControl > 1 && numTreat > 1) {
  pdf("diffexpr-log2correlation.pdf")
  par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1.75,cex.axis=1.75)
  plotFoldChangeCorrelation(resdata,axisLim=4,numControl=numControl)
  dev.off()
}

# heatmap
df <- resdata
# df <- df[df$padj < 0.05, ]
df <- na.omit(df)
dataTable <- as.matrix(df[ ,13:dim(df)[2]])
head(dataTable)
dataTable <- log2(dataTable / rowMeans(dataTable))
dataTable[!is.finite(dataTable)] <- NA
dataTable <- na.omit(dataTable)

clust <- kmeans(dataTable,3)
dataTable <- dataTable[order(clust$cluster), ]
labels <- df$label[order(clust$cluster)]

# colnames(dataTable) <- gsub("_reads","",colnames(dataTable))

pdf("diffexpr-log2heatmap.pdf",useDingbats=F,height=5,width=4)
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=0.75,cex.axis=0.75)
my_palette <- colorRampPalette(c("#0071bc", "gray95", "#d85218"))(n = 255)
heatmap.2(dataTable, col=my_palette, margins=c(7,5), labRow=NA, cexCol=0.75, srtCol=45,
          density.info="none", trace="none", key=T, key.title="NA", key.xlab="log2 over mean",
          dendrogram="none", symkey=T, scale="none", Rowv=F, Colv=F, breaks=seq(-1,1,length.out=256))
#key.xtickfun=function() { return(list(at=c(-0,0.5,1),labels=as.character(c(-1,0,1)))) })
title(main=list(paste0(dim(dataTable)[1]," selected genes"), cex=0.75))

dev.off()

