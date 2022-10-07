#!/usr/bin/env Rscript

library(hodgeslabR)
library(pheatmap)
library(reshape)
library(Seurat)

numSamples <- 2
laplaceSmooth <- 0.01

# run after 1_preprocess_Seurat.R
load("All_conditions_SCTransform.Rda")

# data is in gbm@assays[["SCT"]]@data

# hch note: all cells are assigned (confirmed):
#   nrow(barcodeIdents) == ncol(gbm@assays[["SCT"]]@data) == TRUE
# hch note: barcode order is also the same (confirmed):
#   sum(barcodeIdents$barcode == colnames(gbm@assays[["SCT"]]@data)) == 21947
# no need to take intersection or reorder to align barcodes

# assign cell types to clusters based on manual curation
identFile <- "cell_clusters.csv"
barcodeIdents <- read.table(identFile, sep=",", header=T)
barcodeIdents$dataSet <- substr(barcodeIdents$barcode,18,18)
Idents(gbm) <- barcodeIdents$cellType

control <- 1
treat <- 2

for(cellType in unique(barcodeIdents$cellType)) {

  expr_c <- rowMeans(as.matrix(gbm@assays[["SCT"]]@data[,barcodeIdents$dataSet == control & barcodeIdents$cellType==cellType]))
  expr_t <- rowMeans(as.matrix(gbm@assays[["SCT"]]@data[,barcodeIdents$dataSet == treat & barcodeIdents$cellType==cellType]))
  
  dfp <- data.frame(control=expr_c, treat=expr_t)
  dfp <- dfp[expr_c > 0 & expr_t > 0, ] # drop undetected genes
  dfp$baseMean <- log10(rowMeans(dfp))
  dfp$log2FoldChange <- log2((dfp$treat+laplaceSmooth)/(dfp$control+laplaceSmooth))
  dfp$label <- row.names(dfp)
  xmin <- as.numeric(quantile(dfp$baseMean,0.10))
  xmax <- max(dfp$baseMean)
  #plot(dfp$baseMean,dfp$log2FoldChange,pch=".",xlim=c(xmin,xmax),ylim=c(-4,4))
  #plot(log2(dfp$control),log2(dfp$treat),pch=".",xlim=c(xmin,xmax),ylim=c(xmin,xmax))
  dfp$size = 0.25
  dfp$color <- "#a0a0a0"
  dfp$alpha <- 1
  sel <- abs(dfp$log2FoldChange) > 1
  dfp$size[sel] <- 0.25
  dfp$color[sel] <- "#000000"
  dfp$alpha[sel] <- 1
  
  # remove mitochondrial genes
  sel <- !grepl("^mt-",row.names(dfp))
  dfp <- dfp[sel,]
  
  dfp$label[dfp$color != "#000000"] <- ""
  p1 <- ggplot(dfp,aes(x=baseMean,y=log2FoldChange,size=size,color=color,alpha=alpha,label=label)) +
    theme_hodgeslab_basic() +
    geom_point(shape=16) +
    geom_hline(yintercept = 0, size=0.25*linescale) +
    geom_hline(yintercept = 1, size=0.25*linescale, linetype=2) +
    geom_hline(yintercept = -1, size=0.25*linescale, linetype=2) +
    scale_alpha_identity() +
    scale_size_identity() +
    scale_color_identity() +
    geom_text(colour = "red") +
    xlim(c(xmin,xmax)) +
    ylim(-3,3) +
    xlab("Log10 mean expression") +
    ylab("Log2 fold change")
  
  #geom_text_repel(size=2, min.segment.length = 0, segment.size = 0.25*linescale)
  
  fn <- paste0("ma_",control,"-",treat,"_",cellType,"_labeled.pdf")
  pdf(fn,height=1.5, width=1.5, useDingbats = F)
  print(p1)
  dev.off()

  fn2 <- paste0("ma_",control,"-",treat,"_",cellType,".txt")
  write.table(dfp, fn2, quote=F, sep="\t", col.names=T, row.names = T)
  
  outSub <- paste0(control,"-",treat,"/",cellType)
  fn3 <- paste0(outSub,"/gseaData.rnk")
  # dfgsea <- data.frame(gene=toupper(dfp$label), log2FoldChange=as.character(dfp$log2FoldChange))
  dfgsea <- data.frame(gene=toupper(row.names(dfp)), log2FoldChange=as.character(dfp$log2FoldChange))
  #  dfgsea <- dfgsea[order(dfgsea$log2FoldChange, decreasing = T), ]
  dfgseah <- data.frame(gene="#NAME", log2FoldChange="log2FoldChange")
  dfgsea <- rbind(dfgseah,dfgsea)
  dir.create(outSub, recursive = T)
  write.table(dfgsea, fn3, quote=F, sep="\t", col.names=F, row.names = F)
}


