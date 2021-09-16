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

# remove mostly empty genes, suggested from Michael Love
# removes the genes which have all 0's except a single sample with a high count
keep <- rowSums(counts(dds) >= 1) >= 3
dds <- dds[keep,]

# do everything automagically
dds <- DESeq(dds, fitType="local")


# or manally
# because ChIP seq has background peak calls, we normalize not to housekeeping genes:
# dds <- estimateSizeFactors(dds)
# dds <- estimateSizeFactors(dds, controlGenes = (rowMeans(countdata) > 1000))
# but instead set the sizeFactors to be the total number of reads within peaks
# sizeFactors(dds) <- colSums(countdata) / 1e6
# sizeFactors(dds) <- colSums(counts(dds)) / 1e6

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

#dds <- estimateDispersions(dds, fitType = "local")
#dds <- nbinomWaldTest(dds)

# plot dispersions
pdf("qc-dispersions.pdf")
par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=1.75,cex.axis=1.75)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

rld <- rlogTransformation(dds, blind=FALSE)
rlddata <- as.data.frame(assay(rld))
# row.names(rlddata) <- c(1:nrow(rlddata))

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

# res <- results(dds)
res <- lfcShrink(dds, coef=2)
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
resdata_rld <- resdata
head(resdata)

## Merge with normalized count data
resdata <- merge(resdata, as.data.frame(counts(dds, normalized=T)), by="row.names", sort=F)
#resdata <- merge(as.data.frame(res), as.data.frame(rlddata), by="row.names", sort=F)
names(resdata)[1] <- "id"
head(resdata)
write.table(resdata, file="diffexpr-results.txt", sep="\t", quote=FALSE, row.names=F)

# output variance stabilizing rlog transformed counts as separate file
resdata_rld <- merge(resdata_rld, as.data.frame(rlddata), by="row.names", sort=F)
names(resdata_rld)[1] <- "id"
write.table(resdata_rld, file="diffexpr-results_rld.txt", sep="\t", quote=FALSE, row.names=F)

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
#DESeq2::plotMA(res, ylim=c(-2.25,2.25), ylab="log2 fold change", xlab="mean counts (norm.)")
maplot(res, ylim=c(-2.25,2.25), ylab="log2 fold change", xlab="log10 counts (norm.)")
abline(h=-log2(1.5),col="#0071bc",lwd=2)
abline(h=log2(1.5),col="#d85218",lwd=2)
#require(calibrate)
#with(subset(res, is.element(label,c("Ezh2","Suz12","Rbbp4","Eed","Rybp","Cbx7","Rnf2","Gapdh","Hsp90","Smarca2","Actb","Smarca4"))),
#  text(log10(baseMean), log2FoldChange, labels=label, cex=1.75, col="#333333", pos=4, offset=0.5))
#with(subset(res, is.element(label,c("Ezh2","Suz12","Rbbp4","Eed","Rybp","Cbx7","Rnf2","Gapdh","Hsp90","Smarca2","Actb","Smarca4"))),
#  points(log10(baseMean), log2FoldChange, col="#aaaaaa", pch=16, cex=0.5))

#axis(side=1,at=c(10,100,1000))
dev.off()

#llplot <- function (res, thresh=0.10, foldthresh=log2(1.5), labelthresh=1.5, labelsig=F, minCount=0, ...) {
#  with(subset(res, (!(padj<thresh & abs(log2FoldChange) > foldthresh) | is.na(padj)) & baseMean > minCount), heatscatter(log10((c1+c2)/2), log10((t1+t2)/2), pch=16, cex=0.5, colpal="greys", rev=F, main="", ...))
#  with(subset(res, padj<thresh & log2FoldChange < -foldthresh & baseMean > minCount), heatscatterpoints(log10((c1+c2)/2), log10((t1+t2)/2), colpal="blues", rev=F, pch=16, cex=0.5))
#  with(subset(res, padj<thresh & log2FoldChange > foldthresh & baseMean > minCount), heatscatterpoints(log10((c1+c2)/2), log10((t1+t2)/2), colpal="oranges", rev=F, pch=16, cex=0.5))
#  if (labelsig) {
#    require(calibrate)
#    with(subset(res, padj<thresh & abs(log2FoldChange) > labelthresh), textxy(baseMean, log2FoldChange, labs=label, cex=1, col="#333333"))
#  }
#}

#llplotc <- function (res, thresh=0.10, foldthresh=log2(1.5), labelthresh=1.5, labelsig=F, minCount=0, ...) {
#  with(subset(res, ( !is.na(padj) | is.na(padj) )), heatscatter(log10(c1), log10(c2), pch=16, cex=0.5, colpal="greys", rev=F, main="", ...))
#}

#llplott <- function (res, thresh=0.10, foldthresh=log2(1.5), labelthresh=1.5, labelsig=F, minCount=0, ...) {
#  with(subset(res, ( !is.na(padj) | is.na(padj) )), heatscatter(log10(t1), log10(t2), pch=16, cex=0.5, colpal="greys", rev=F, main="", ...))
#}

resll <- resdata
colnames(resll)[c(13,14,13+numControl,13+numControl+1)] <- c("c1","c2","t1","t2")

#pdf("diffexpr-llplot.pdf",useDingbats=F)
#par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=2.25,cex.axis=2.25)
#llplot(resll, ylab="log10 treatment expr. (norm.)", xlab="log10 control expr. (norm.)")
#with(subset(res, is.element(label,c("Ezh2","Suz12","Rbbp4","Eed","Rybp","Cbx7","Rnf2","Gapdh","Hsp90","Smarca2","Actb","Smarca4"))),
#  text(log10(baseMean), log2FoldChange, labels=label, cex=1.75, col="#333333", pos=4, offset=0.5))
#with(subset(res, is.element(label,c("Ezh2","Suz12","Rbbp4","Eed","Rybp","Cbx7","Rnf2","Gapdh","Hsp90","Smarca2","Actb","Smarca4"))),
#  points(log10(baseMean), log2FoldChange, col="#aaaaaa", pch=16, cex=0.5))
#dev.off()


#pdf("diffexpr-llplot_control.pdf",useDingbats=F)
#par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=2.25,cex.axis=2.25)
#llplotc(resll, ylab="log10 control rep2 (norm.)", xlab="log10 control rep1 (norm.)")
#with(subset(res, is.element(label,c("Ezh2","Suz12","Rbbp4","Eed","Rybp","Cbx7","Rnf2","Gapdh","Hsp90","Smarca2","Actb","Smarca4"))),
#  text(log10(baseMean), log2FoldChange, labels=label, cex=1.75, col="#333333", pos=4, offset=0.5))
#with(subset(res, is.element(label,c("Ezh2","Suz12","Rbbp4","Eed","Rybp","Cbx7","Rnf2","Gapdh","Hsp90","Smarca2","Actb","Smarca4"))),
#  points(log10(baseMean), log2FoldChange, col="#aaaaaa", pch=16, cex=0.5))
#dev.off()


#pdf("diffexpr-llplot_treat.pdf",useDingbats=F)
#par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=2.25,cex.axis=2.25)
#llplott(resll, ylab="log10 treatment rep2 (norm.)", xlab="log10 treatment rep1 (norm.)")
#with(subset(res, is.element(label,c("Ezh2","Suz12","Rbbp4","Eed","Rybp","Cbx7","Rnf2","Gapdh","Hsp90","Smarca2","Actb","Smarca4"))),
#  text(log10(baseMean), log2FoldChange, labels=label, cex=1.75, col="#333333", pos=4, offset=0.5))
#with(subset(res, is.element(label,c("Ezh2","Suz12","Rbbp4","Eed","Rybp","Cbx7","Rnf2","Gapdh","Hsp90","Smarca2","Actb","Smarca4"))),
#  points(log10(baseMean), log2FoldChange, col="#aaaaaa", pch=16, cex=0.5))
#dev.off()



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

# heatmap
# df <- resdata_rld
# df <- df[df$padj < 0.25, ]
# df <- na.omit(df)
# dataTable <- as.matrix(df[ ,13:dim(df)[2]])
# head(dataTable)
# #dataTable <- log2(dataTable / rowMeans(dataTable))
# dataTable <- t(scale(t(dataTable)))
# dataTable[!is.finite(dataTable)] <- NA
# dataTable <- na.omit(dataTable)

# clust <- kmeans(dataTable,2)
# dataTable <- dataTable[order(clust$cluster), ]
# labels <- df$label[order(clust$cluster)]

# colnames(dataTable) <- gsub("_reads","",colnames(dataTable))

# pdf("diffexpr-log2heatmap.pdf",useDingbats=F,height=5,width=4)
# par(oma=c(1,1,1,1),mar=c(5,5,4,1)+0.1,cex.lab=0.75,cex.axis=0.75)
# my_palette <- colorRampPalette(c("#0071bc", "gray95", "#d85218"))(n = 255)
# heatmap.2(dataTable, col=my_palette, margins=c(7,5), labRow=labels, cexCol=0.75, srtCol=45,
#           density.info="none", trace="none", key=T, key.title="NA", key.xlab="Z score",
#           dendrogram="none", symkey=T, scale="none", Rowv=F, Colv=F, breaks=seq(-1,1,length.out=256))
# #key.xtickfun=function() { return(list(at=c(-0,0.5,1),labels=as.character(c(-1,0,1)))) })
# title(main=list(paste0(dim(dataTable)[1]," selected genes"), cex=0.75))

# dev.off()

