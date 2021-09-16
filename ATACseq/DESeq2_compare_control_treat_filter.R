#!/usr/bin/env Rscript

# assumes will read in tab-delimited file with one header line, with columns as the following:
# id
# chr
# start
# end
# label
# control1_counts
# control2_counts
# ...
# treat1_counts
# treat2_counts
# ...

# command line args
args <- commandArgs(TRUE)
tableFile <- args[1]
numControl <- as.numeric(args[2])
# numTreat <- as.numeric(args[3])

# organize tabulated data
dataTable  <- read.table(tableFile, header=TRUE, row.names=1, sep="\t")
dataChr    <- dataTable[ ,1]
dataStart  <- dataTable[ ,2]
dataEnd    <- dataTable[ ,3]
dataLabel  <- dataTable[ ,4]
dataCounts <- dataTable[ ,5:ncol(dataTable)]
# colnames(dataCounts) <- gsub("\\_reads$", "", colnames(dataCounts))
countdata  <- as.matrix(dataCounts)
countdata[countdata < 0] <- 0
head(countdata)

numTreat <- ncol(countdata) - numControl

# set conditions for each dataset
condition <- factor(c(rep("control", numControl), rep("treat", numTreat)))
#replicate <- factor(c(seq(1,numControl,1), seq(1,numTreat,1)))

# load DESeq2 and process data
library(DESeq2)
coldata <- data.frame(row.names=colnames(countdata), condition)
coldata$dataset <- colnames(countdata)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds$condition <- relevel(dds$condition, "control")

# remove mostly empty genes, suggested from Michael Love
# removes the genes which have all mostly high counts but single
# samples with  0's
keep <- rowSums(counts(dds) >= 1) >= 3
dds <- dds[keep,]
dataLabel <- dataLabel[keep]
dataChr <- dataChr[keep]
dataStart <- dataStart[keep]
dataEnd <- dataEnd[keep]

# do everything automagically
# dds <- DESeq(dds, fitType="local")
# or manally
# because ChIP seq has background peak calls, we normalize not to housekeeping genes:
dds <- estimateSizeFactors(dds)
# dds <- estimateSizeFactors(dds, controlGenes = (rowMeans(countdata) > 1000))
# but instead set the sizeFactors to be the total number of reads within peaks
# sizeFactors(dds) <- colSums(countdata) / 1e6

# output sizeFactors for later analyses
write.table(data.frame(value=sizeFactors(dds)),"sizeFactors.txt",sep="\t",quote=F)

dds <- estimateDispersions(dds, fitType = "local")
dds <- nbinomWaldTest(dds)

# plot dispersions
pdf("qc-dispersions.pdf",useDingbats=F)
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1.75,cex.axis=1.75)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

rld <- rlogTransformation(dds, blind=FALSE)
rlddata <- as.data.frame(assay(rld), row.names=c(1:nrow(countdata)))
head(assay(rld))
# hist(assay(rld))
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# plot heatmap distance of samples
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
pdf("qc-heatmap-samples.pdf",useDingbats=F)
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1,cex.axis=1)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margins=c(12, 12), main="Sample Distance Matrix", cexRow=1,cexCol=1,srtCol=45)
dev.off()

# plot PCA
pdf("qc-pca.pdf",useDingbats=F,width=12,height=7)
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1.75,cex.axis=1.75)
DESeq2::plotPCA(rld, intgroup=c("condition", "dataset"))
dev.off()

# res <- results(dds)
res <- lfcShrink(dds, coef=2)
res$label <- dataLabel
res$chr <- dataChr
res$startPos <- dataStart
res$endPos <- dataEnd
table(res$padj<0.05)
## Order by adjusted p-value
# res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdatarld <- merge(as.data.frame(res), as.data.frame(rlddata), by="row.names", sort=FALSE)
names(resdata)[1] <- "id"
head(resdata)

write.table(resdata, file="diffexpr-results.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(resdatarld, file="diffexpr-results_rld.txt", sep="\t", quote=FALSE, row.names=FALSE)

## output MA plot
#maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
#  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
#  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
#  if (labelsig) {
#    require(calibrate)
#    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=id, cex=textcx, col=2))
#  }
#}

library(LSD)
maplot <- function (res, thresh=0.10, foldthresh=log2(1.5), labelthresh=1.5, labelsig=F, minCount=0, ...) {

#  with(subset(res, (!(padj<thresh & abs(log2FoldChange) > foldthresh) | is.na(padj)) & baseMean > minCount),
#    heatscatter(log10(baseMean), log2FoldChange, pch=16, cex=0.5, colpal="greys", rev=F, main="", ...))

  sel <- (!(res$padj<thresh & abs(res$log2FoldChange) > foldthresh) | is.na(res$padj)) & res$baseMean > minCount
  if (sum(sel, na.rm=T) > 0)
    heatscatter(log10(res$baseMean[sel]), res$log2FoldChange[sel], pch=16, cex=0.5, colpal="greys", rev=F, main="", ...)

  sel <- !is.na(res$padj) & res$padj<thresh & res$log2FoldChange < -foldthresh & res$baseMean > minCount
  if (sum(sel, na.rm=T) == 1)
    points(log10(res$baseMean[sel]), res$log2FoldChange[sel], col=RColorBrewer::brewer.pal(3,"Blues")[3], pch=16, cex=0.5)
  else if (sum(sel, na.rm=T) > 0)
    heatscatterpoints(log10(res$baseMean[sel]), res$log2FoldChange[sel], colpal="blues", rev=F, pch=16, cex=0.5)

  sel <- !is.na(res$padj) & res$padj<thresh & res$log2FoldChange > foldthresh & res$baseMean > minCount
  if (sum(sel, na.rm=T) == 1)
    points(log10(res$baseMean[sel]), res$log2FoldChange[sel], col=RColorBrewer::brewer.pal(3,"Oranges")[3], pch=16, cex=0.5)
  else if (sum(sel, na.rm=T) > 1)
    heatscatterpoints(log10(res$baseMean[sel]), res$log2FoldChange[sel], colpal="oranges", rev=F, pch=16, cex=0.5)

  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh & abs(log2FoldChange) > labelthresh),
      textxy(baseMean, log2FoldChange, labs=label, cex=1, col="#333333"))
  }
}


pdf("diffexpr-maplot.pdf",useDingbats=F)
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=2.25,cex.axis=2.25)
#DESeq2::plotMA(res, ylim=c(-2,2), ylab="log2 fold change", xlab="mean counts (norm.)")
maplot(res, ylim=c(-2.25,2.25), ylab="log2 fold change", xlab="log10 peak intensity (norm.)")
abline(h=-log2(1.5),col="#0071bc",lwd=2)
abline(h=log2(1.5),col="#d85218",lwd=2)
#require(calibrate)
#with(subset(res, is.element(label,c("Ezh2","Rybp","Suz12","Gapdh","Hsp90","Smarca2","Actb","Rnf2","Smarca4"))),
#  text(log10(baseMean), log2FoldChange, labels=label, cex=1.75, col="#333333", pos=4, offset=0.5))

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
pdf("diffexpr-volcanoplot.pdf",useDingbats=F,width=12,height=7)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8) # , xlim=c(-5, 5), ylim=c(-1, 55))
dev.off()

library(boot)

corbootf <- function(d,i) {
	d2 <- d[i, ]
	return(cor(d2$rep1, d2$rep2, method="pearson"))
}

plotFoldChangeCorrelation <- function (resdata, axisLim=1,sigthresh=0.10, numControl=2, numTreat=2, main=NULL, logthresh=log2(1.5)) {
	df <- as.data.frame(resdata[ ,c(7,12,13,12+numControl,13+numControl)])
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

	sigPoints <- df$padj < sigthresh & abs(resdata$log2FoldChange) > logthresh
	sigPointsIncr <- !is.na(df$padj) & df$padj < sigthresh & resdata$log2FoldChange > logthresh
	sigPointsDecr <- !is.na(df$padj) & df$padj < sigthresh & resdata$log2FoldChange < -logthresh
	nsigPoints <- ! sigPoints | is.na(df$padj)
	points(df$rep1[nsigPoints], df$rep2[nsigPoints], pch=16, col=rgb(0,0,0,0.2),cex=0.5)

	if (sum(sigPointsIncr, na.rm=T) == 1)
	  points(df$rep1[sigPointsIncr], df$rep2[sigPointsIncr], col=RColorBrewer::brewer.pal(3,"Oranges")[3], pch=16, cex=0.5)
	else if(sum(sigPointsIncr, na.rm=T) > 1)
	  heatscatterpoints(df$rep1[sigPointsIncr], df$rep2[sigPointsIncr], colpal="oranges", rev=F, pch=16, cex=0.5)

	if (sum(sigPointsDecr, na.rm=T) == 1)
	  points(df$rep1[sigPointsDecr], df$rep2[sigPointsDecr], col=RColorBrewer::brewer.pal(3,"Blues")[3], pch=16, cex=0.5)
	else if(sum(sigPointsDecr, na.rm=T) > 1)
	  heatscatterpoints(df$rep1[sigPointsDecr], df$rep2[sigPointsDecr], colpal="blues", rev=F, pch=16, cex=0.5)

	rval = cor(df$rep1, df$rep2, method="pearson")
	
	# bootstrapping with 100 replications 
	results <- boot(data=df, statistic=corbootf, R=100)
	rsqci <- boot.ci(results, type="basic")
	label <- sprintf("R = %.3f\n(%.3f, %.3f) 95%% CI",rval,rsqci$basic[4],rsqci$basic[5],3)
	text(-0.5*axisLim,0.75*axisLim,label,offset=0,cex=1.4)
}

if(numControl > 1 && numTreat > 1) {
  pdf("diffexpr-log2correlation.pdf", useDingbats=F)
  par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=2.25,cex.axis=2.25)
  plotFoldChangeCorrelation(resdata,axisLim=2.25,numControl=numControl)
  dev.off()
}
