# script to perform pseudo-bulk DGA
# setwd("~/Desktop/DS_training_session_031324/scripts_and_datasets")

library(ExperimentHub)
library(Seurat)
library(DESeq2)
library(tidyverse)


# get data
eh <- ExperimentHub()
query(eh, "Kang")

sce <- eh[["EH2259"]]
seu.obj <- as.Seurat(sce, data = NULL)
View(seu.obj@meta.data)

# QC and filtering
# explore QC

# get mito percent
seu.obj$mitoPercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')
View(seu.obj@meta.data)

# filter
seu.filtered <- subset(seu.obj, subset = nFeature_originalexp > 200 & nFeature_originalexp < 2500 &
         nCount_originalexp > 800 & 
         mitoPercent < 5 &
         multiplets == 'singlet')

# run Seurat's standard workflow steps
seu.filtered <- NormalizeData(seu.filtered)
seu.filtered <- FindVariableFeatures(seu.filtered)
seu.filtered <- ScaleData(seu.filtered)
seu.filtered <- RunPCA(seu.filtered)
ElbowPlot(seu.filtered)
seu.filtered <- RunUMAP(seu.filtered, dims = 1:20)

# visualize 
cell_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'cell', label = TRUE)
cond_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'stim')

cell_plot|cond_plot

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level


# convert to data.frame

# get values where to split

# split data.frame

# fix colnames and transpose



# Let's run DE analysis with B cells
# 1. Get counts matrix

# 2. generate sample level metadata

# get more information from metadata


# perform DESeq2 --------
# Create DESeq2 object

# filter

# run DESeq2

# Check the coefficients for the comparison

# Generate results object


