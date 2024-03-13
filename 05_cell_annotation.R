# script to identify cluster identity
# Finding markers in every cluster
# Finding conserved markers (optional)
# Finding markers DE between conditions
# setwd("~/Desktop/DS_training_session_031324/scripts_and_datasets")
set.seed(42)

library(Seurat)
library(tidyverse)
#install.packages('BiocManager')
#BiocManager::install('multtest')
#BiocManager::install("DESeq2")
#install.packages('metap')

# load data
ifnb_harmony <- readRDS('./datasets/ifnb_harmony.rds')
str(ifnb_harmony)
View(ifnb_harmony@meta.data)

# visualize data

# findAll markers -----------------

# findConserved markers -------------
# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# DefaultAssay(ifnb_harmony) <- 'RNA'


# let's visualize top features



# rename cluster 3 ident


# cells already have annotations provided in the metadata


# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)



# findMarkers between conditions ---------------------

# find markers

# plotting conserved features vs DE features between conditions


