#!/usr/bin/env Rscript

library(Seurat)
library(DoubletFinder) # devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(cowplot)

expFreq <- 0.10 # expect 25% doublets

# multiome data matrix has multiple modalities, dat is a list of matrices for this genome
# names are dat$`Gene Expression` and dat$Peaks

dat <- Seurat::Read10X_h5(filename = "outs/filtered_feature_bc_matrix.h5", use.names = T)

# mitochondrial genes
# automatically
gex <- CreateSeuratObject(counts = dat, assay = "RNA")
gex <- UpdateSeuratObject(gex)
gex <- PercentageFeatureSet(gex, "^mt-", col.name = "percent_mito")
# manually
# total_counts_per_cell <- colSums(dat@assays$RNA@counts)
# mito_genes <- rownames(dat)[grep("^mt-", rownames(dat))]
# dat$percent_mito <- colSums(dat@assays$RNA@counts[mito_genes, ])/total_counts_per_cell

# ribosomal genes
# automatically
gex <- PercentageFeatureSet(gex, "^Rp[sl]", col.name = "percent_ribo")
# manual
# ribo_genes <- rownames(dat)[grep("^Rp[sl]", rownames(dat))]
# dat$percent_ribo <- colSums(dat@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell

# plot QC metrics

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
pdf("qc_DoubletFinder_features.pdf", height=8, width=8, useDingbats=FALSE)
VlnPlot(gex, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) +
    NoLegend()
dev.off()

pdf("qc_DoubletFinder_features_count_vs_genes.pdf", height=4, width=8, useDingbats=FALSE)
FeatureScatter(gex, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
dev.off()

# filter based on percent mito and ribo reads
# note these can be arbitrary thresholds
selected_mito <- WhichCells(gex, expression = percent_mito < 15)
selected_ribo <- WhichCells(gex, expression = percent_ribo > 0.1)

# subset to keep only those cells
gex.filt <- subset(gex, cells = selected_mito)
gex.filt <- subset(gex.filt, cells = selected_ribo)

# plot post QC
pdf("qc_DoubletFinder_features_post.pdf", height=8, width=8, useDingbats=FALSE)
VlnPlot(gex.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) +
    NoLegend()
dev.off()

# run DoubletFinder

# preprocess
gex.filt <- NormalizeData(gex.filt, verbose = F)
gex.filt <- FindVariableFeatures(gex.filt, verbose = F)
gex.filt <- ScaleData(gex.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"), verbose = F)
gex.filt <- RunPCA(gex.filt, verbose = F, npcs = 20)
gex.filt <- RunUMAP(gex.filt, dims = 1:10, verbose = F)

nExp <- round(ncol(gex.filt) * expFreq)
gex.filt <- doubletFinder_v3(gex.filt, pN = 0.25, pK = 0.01, nExp = nExp, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(gex.filt@meta.data)[grepl("DF.classification", colnames(gex.filt@meta.data))]

pdf("qc_DoubletFinder_doublet_detection_umap.pdf", height=4, width=8, useDingbats=FALSE)
cowplot::plot_grid(ncol = 2, DimPlot(gex.filt, group.by = "orig.ident") + NoAxes(),
    DimPlot(gex.filt, group.by = DF.name) + NoAxes())
dev.off()

pdf("qc_DoubletFinder_doublet_detection_genes.pdf", height=4, width=8, useDingbats=FALSE)
VlnPlot(gex.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
dev.off()

# write out names of singlets to file
singlets_gex <- colnames(gex.filt)[gex.filt@meta.data[, DF.name] == "Singlet"]
write.table(singlets_gex, file="qc_DoubletFinder_pass_singlets.txt", quote=F, row.names=F, col.names = F)





