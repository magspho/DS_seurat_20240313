# script to annotate cell types from 20k Human PBMCs from a healthy female donor
# setwd("~/Desktop/DS_training_session_031324/scripts_and_datasets")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

# Input Data 10X CellRanger .HDF5 format --------------
hdf5_obj <- Read10X_h5(filename = './datasets/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)

# QC and Filtering -----------
# explore QC
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)


# It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.


# pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# running steps above to get clusters
View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')

# get reference data -----------

# expression values are log counts (log normalized counts)


# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

# Annotation diagnostics ----------

# ...Based on the scores within cells -----------


# ...Based on deltas across cells ----------


# ...Comparing to unsupervised clustering ------------


