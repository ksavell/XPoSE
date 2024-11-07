# 
library(anndata)
library(tidyverse)
library(Seurat)


library(reticulate)

# only do this part once!
#version <- "3.9.12"
#install_python(version)
#virtualenv_create("my-environment", version = version)
#use_virtualenv("my-environment")
#py_install("anndata", envname = "my-environment", method = "virtualenv")

anndata <- import("anndata")

load("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_September2024/DataForFigures/Robjects/combined08032023.RData")

# pull counts from object
obj <- combined

obj_counts <- (obj[["RNA"]]$counts)
obj_counts <- t(as.matrix(obj_counts))

counts <- obj_counts
genes   <- colnames(counts)
samples <- rownames(counts)
sparse_counts <- as(counts, "dgCMatrix")  # Convert to sparse matrix, if not already

countAD <- AnnData(X   = sparse_counts,   # Create the anndata object
                   var = data.frame(genes=genes,row.names=genes),
                   obs = data.frame(samples=samples,row.names=samples))
write_h5ad(countAD, "combinedcounts.h5ad") # Write it out as h5ad
