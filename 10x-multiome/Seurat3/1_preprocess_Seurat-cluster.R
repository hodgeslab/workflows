#!/usr/bin/env Rscript

# https://satijalab.org/seurat/v3.2/sctransform_vignette.html

# rm(list=ls())

library(Seurat)
library(sctransform)

dat <- Seurat::Read10X_h5(filename = "outs/filtered_feature_bc_matrix.h5", use.names = T)

# mitochondrial genes
# automatically
gbm <- CreateSeuratObject(counts = dat$`Gene Expression`, assay = "RNA")

# store mitochondrial percentage in object meta data
gbm <- PercentageFeatureSet(gbm, pattern = "^mt-", col.name = "percent.mt")

# run sctransform
gbm <- SCTransform(gbm, vars.to.regress = "percent.mt", verbose = FALSE)

# data is in gbm@assays[["SCT"]]@data

save.image("Seurat_SCTransform.Rda")


