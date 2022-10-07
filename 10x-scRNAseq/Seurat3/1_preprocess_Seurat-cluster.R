#!/usr/bin/env Rscript

# https://satijalab.org/seurat/v3.2/sctransform_vignette.html

# rm(list=ls())

library(Seurat)
library(sctransform)

pipestancePath <- "../All_conditions_singlets/outs/filtered_feature_bc_matrix"

gbmData <- Read10X(pipestancePath)
gbm <- CreateSeuratObject(counts = gbmData)

# store mitochondrial percentage in object meta data
gbm <- PercentageFeatureSet(gbm, pattern = "^mt-", col.name = "percent.mt")

# run sctransform
gbm <- SCTransform(gbm, vars.to.regress = "percent.mt", verbose = FALSE)

# data is in gbm@assays[["SCT"]]@data

save.image("All_conditions_SCTransform.Rda")


