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
# processed seurat object -- could also just be a sparse matrix
rds_obj <- readRDS('datasets/ifnb_harmony.rds')
str(rds_obj)

# .loom files
loom_obj <- Connect(filename = "datasets/adult-hem-organs-10X-bone-marrow.loom", mode = 'r')
seurat_loom <- as.Seurat(loom_obj)

# .h5ad format 
# step 1: convert AnnData object to an h5Seurat file
Convert("datasets/adata_SS2_for_download.h5ad", dest = "h5seurat", overwrite = TRUE)

# step 2: Load h5Seurat file into a Seurat object 
seurat_anndata <- LoadH5Seurat("datasets/adata_SS2_for_download.h5seurat")

# .mtx file
mtx_obj <- ReadMtx(mtx = "datasets/raw_feature_bc_matrix/matrix.mtx.gz",
                   features = "datasets/raw_feature_bc_matrix/features.tsv.gz",
                   cells = "datasets/raw_feature_bc_matrix/barcodes.tsv.gz")
pbmc.seurat <- CreateSeuratObject(counts = mtx_obj)

# 10X CellRanger .HDF5 format 
nsclc.sparse.m <- Read10X_h5(filename = "datasets/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
str(nsclc.sparse.m)
cts <- nsclc.sparse.m$`Gene Expression`
cts[1:10,1:10]

# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 29552 features across 42081 samples
