# script to demonstrate reading single cell matrices in various format 
# and converting to seurat object 
# setwd("~/Desktop/DS_training_session_031324/scripts_and_datasets")

# load libraries
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(SeuratDisk)

# .RDS format

# .mtx file

# .loom files

# .h5ad format 
# step 1: convert AnnData object to an h5Seurat file

# step 2: Load h5Seurat file into a Seurat object 


# 10X CellRanger .HDF5 format 


# Initialize the Seurat object with the raw (non-normalized data).
