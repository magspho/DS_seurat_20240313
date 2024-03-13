# script to integrate across conditions using Harmony

# set seed for reproducibility
set.seed(42)

library(harmony)
library(Seurat)
# install.packages("devtools")
# library(devtools)
# devtools::install_github('satijalab/seurat-data',force = TRUE)
library(SeuratData)
library(tidyverse)
library(ggplot2)

# get data -------------------------
AvailableData()
# install dataset
InstallData("ifnb")
# load dataset
LoadData("ifnb")
str(ifnb)
# QC and filtering
ifnb[["percent.mt"]] <- PercentageFeatureSet(object=ifnb, pattern = '^MT-')
View(ifnb@meta.data)
# explore QC

# filter
ifnb
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          percent.mt < 5)

# standard workflow steps
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
all.genes <- rownames(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered, features=all.genes)
ifnb.filtered <- RunPCA(ifnb.filtered, features=VariableFeatures(object = ifnb.filtered))
ElbowPlot(ifnb.filtered)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')


# run Harmony -----------

# Do UMAP and clustering using ** Harmony embeddings instead of PCA **

# visualize 
