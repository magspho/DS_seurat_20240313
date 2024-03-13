# filter out doublets: DoubletFinder
# setwd("~/Desktop/DS_training_session_031324/scripts_and_datasets")

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

# create counts matrix
#cts <- ReadMtx(mtx = 'datasets/raw_feature_bc_matrix/matrix.mtx.gz',
#        features = 'datasets/raw_feature_bc_matrix/features.tsv.gz',
#        cells = 'datasets/raw_feature_bc_matrix/barcodes.tsv.gz')
# create Seurat object
#pbmc.seurat <- CreateSeuratObject(counts = cts)

# QC and Filtering
pbmc.seurat <- PercentageFeatureSet(pbmc.seurat, 
                                    pattern = "^MT-", # Uses all gene names that begin with "MT-"
                                    col.name = "percent.mt")        
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 percent.mt < 10)

pbmc.seurat
pbmc.seurat.filtered


# pre-process standard workflow
pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(pbmc.seurat.filtered, npcs=100)
ElbowPlot(pbmc.seurat.filtered, ndims = 100)
pbmc.seurat.filtered <- FindNeighbors(pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(pbmc.seurat.filtered, dims = 1:20)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_pbmc <- paramSweep(pbmc.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
pbmc.seurat.filtered <- doubletFinder(pbmc.seurat.filtered, 
                                      PCs = 1:20, 
                                      pN = 0.25, 
                                      pK = pK, 
                                      nExp = nExp_poi.adj,
                                      reuse.pANN = FALSE, sct = FALSE)
View(pbmc.seurat.filtered@meta.data)
# visualize doublets
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.26_691")

# number of singlets and doublets
table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_0.26_691)
