# XPoSE-seq: export main and proof-of-concept datasets for MapMyCells
#
# Original MapMyCells export workflow authors: Katherine E. Savell and Padmashri Saravanan
# Adaptation and consolidation for final repository: Katherine E. Savell
#
# Run from the base XPoSE repository directory.
# This script exports one sparse h5ad file per nucleus capture for annotation
# with Allen Institute MapMyCells.
#
# MapMyCells: https://portal.brain-map.org/atlases-and-data/bkp/mapmycells

rm(list = ls())

# First-time Python setup --------------------------------------------------
# Uncomment this block once when setting up a new machine, then comment it
# out again before subsequent runs.
#
# library(reticulate)
# python_version <- "3.9.12"
# install_python(python_version)
# virtualenv_create("my-environment", version = python_version)
# use_virtualenv("my-environment")
# py_install("anndata", envname = "my-environment", method = "virtualenv")

# Load packages -----------------------------------------------------------
library(Seurat)
library(reticulate)

use_virtualenv("my-environment")
anndata <- import("anndata")

# Paths -------------------------------------------------------------------
dataset_files <- c(
  main = "output/01_metadata_clustering_qc/main_metadata_qc.rds",
  poc  = "output/01_metadata_clustering_qc/poc_metadata_qc.rds"
)

output_dir <- "output/01_metadata_clustering_qc/mapmycells"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Export one h5ad per nucleus capture -------------------------------------
for (dataset_name in names(dataset_files)) {
  
  seur_obj <- readRDS(dataset_files[[dataset_name]])
  seur_obj <- JoinLayers(seur_obj, assay = "RNA")
  
  dataset_output_dir <- file.path(output_dir, dataset_name)
  dir.create(dataset_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  captures <- unique(seur_obj$capture)
  
  for (capture_name in captures) {
    
    capture_cells <- colnames(seur_obj)[seur_obj$capture == capture_name]
    capture_obj <- subset(seur_obj, cells = capture_cells)
    
    counts_mat <- LayerData(
      capture_obj,
      assay = "RNA",
      layer = "counts"
    )
    
    genes <- rownames(counts_mat)
    cell_ids <- colnames(counts_mat)
    
    # Keep the expression matrix sparse during transposition/export.
    sparse_counts <- methods::as(t(counts_mat), "dgCMatrix")
    
    count_ad <- anndata$AnnData(
      X = sparse_counts,
      var = data.frame(
        gene = genes,
        row.names = genes
      ),
      obs = data.frame(
        cell_id = cell_ids,
        row.names = cell_ids
      )
    )
    
    count_ad$write_h5ad(
      file.path(dataset_output_dir, paste0(capture_name, ".h5ad"))
    )
    
    rm(capture_obj, counts_mat, sparse_counts, count_ad)
    invisible(gc())
  }
  
  rm(seur_obj)
  invisible(gc())
}
