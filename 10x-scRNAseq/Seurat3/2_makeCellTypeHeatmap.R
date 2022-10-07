#!/usr/bin/env Rscript

library(hodgeslabR)
# library(pheatmap)
library(Seurat)

numSamples <- 2
identFile <- "cell_clusters.csv"

# run after 1_preprocess_Seurat.R
load("All_conditions_SCTransform.Rda")
# data is in gbm@assays[["SCT"]]@data

# hch note: all cells are assigned (confirmed):
#   nrow(barcodeIdents) == ncol(gbm@assays[["SCT"]]@data) == TRUE
# hch note: barcode order is also the same (confirmed):
#   sum(barcodeIdents$barcode == colnames(gbm@assays[["SCT"]]@data)) == 21947
# no need to take intersection or reorder to align barcodes

# assign cell types to clusters based on manual curation #Do this!
barcodeIdents <- read.table(identFile, sep=",", header=T)
barcodeIdents$dataSet <- substr(barcodeIdents$barcode,18,18)

# remove erythroid cells #? TRY BOTH WAYS
# sel <- c("Erythroid")
# barcodeIdents <- barcodeIdents[! barcodeIdents$cellType %in% sel, ]

# select only annotated cells and assign cell type identities
gbm <- gbm[,barcodeIdents$barcode]
Idents(gbm) <- barcodeIdents$cellType

# find markers for every cluster compared to all remaining cells, report only the positive ones
gbm.markers <- FindAllMarkers(gbm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gbm.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# take top 6 most distinguishing features in each cluster for heatmap
# hch note: heatmap is raster by default
topN <- gbm.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_logFC)

# below is rather complex. this is done because some cells are very lowly
# represented. we do this complex procedure to scale the number of cells in
# each cell type to be represented logarithmically
#for(seldataSet in 1:numSamples) {
for(seldataSet in 1:numSamples) {
  # find least frequent cell type in each data set and log2 scale everything according to it
  classSampleSize <- data.frame(cellType=unique(barcodeIdents$cellType),
                                numSelect = sapply(unique(barcodeIdents$cellType),
                                                   function(x) {
                                                     sizeFactor <- sum(barcodeIdents$cellType == x & barcodeIdents$dataSet == seldataSet) }))
  minSample <- min(classSampleSize$numSelect)

  # rescale according to log2 scale
  classSampleSize <- data.frame(cellType=unique(barcodeIdents$cellType),
                                numSelect = sapply(unique(barcodeIdents$cellType),
                                                   function(x) {
                                                     sizeFactor <- sum(barcodeIdents$cellType == x & barcodeIdents$dataSet == seldataSet)
                                                     sizeFactor <- minSample + floor(minSample*log2((sizeFactor/minSample)+1))
                                                     #sizeFactor <- minSample + floor(minSample*log2((sizeFactor/minSample)+1)) #If log2(0) undefined is a problem, try adding +1
                                                   }))
  cellSelect <- sapply(unique(barcodeIdents$cellType),
                       function(x) {
                         selData <- barcodeIdents[barcodeIdents$cellType == x & barcodeIdents$dataSet == seldataSet,]
                         position <- 1:nrow(selData)
                         selBarcodes <- as.character(selData$barcode[1:min(nrow(selData),
                                                                           classSampleSize$numSelect[classSampleSize$cellType == x])])
                       })
  cellSelect <- unlist(cellSelect)
  
  outfn <- paste0("summary_cellType_log2_heatmap_",seldataSet,".pdf")
  pdf(outfn, height=5, width=5, useDingbats = F)
  p1 <- DoHeatmap(gbm, cells=cellSelect, assay="SCT",
            features = topN$gene, raster=F) + #NoLegend() + 
    # default colors:
    #  scale_fill_gradientn(colors = c("magenta", "black", "yellow"))
    # scale_fill_viridis()
    scale_fill_gradientn(colors = c("#8b008b", "black", "yellow"))
  print(p1)
  dev.off()
  
  
  outfn <- paste0("summary_cellType_linear_heatmap_",seldataSet,".pdf")
  pdf(outfn, height=5, width=5, useDingbats = F)
  p1 <- DoHeatmap(gbm, assay="SCT",
                  features = topN$gene, raster=F) + #NoLegend() + 
    # default colors:
    #  scale_fill_gradientn(colors = c("magenta", "black", "yellow"))
    # scale_fill_viridis()
    scale_fill_gradientn(colors = c("#8b008b", "black", "yellow"))
  print(p1)
  dev.off()
}

