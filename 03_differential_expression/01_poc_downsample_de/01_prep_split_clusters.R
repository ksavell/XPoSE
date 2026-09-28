#!/usr/bin/env Rscript
# XPoSE-seq: prepare POC clusters for downsampling differential expression
#
# Run from the base XPoSE repository directory.
# This one-time preparation step splits the final POC combined Seurat object
# into one RDS file per cell type for the Biowulf downsampling workflow.

suppressPackageStartupMessages({
  library(Seurat)
})

# Paths -------------------------------------------------------------------
input_file <- "output/01_metadata_clustering_qc/poc_combined_annotated.rds"
clusters_file <- "03_differential_expression/01_downsample_de/clusters_kept.txt"
data_root <- "output/03_differential_expression/02_downsample_de/hpc_input"

dir.create(data_root, showWarnings = FALSE, recursive = TRUE)

# Load data ---------------------------------------------------------------
seur_obj <- readRDS(input_file)
clusters <- readLines(clusters_file)
clusters <- trimws(clusters)
clusters <- clusters[nzchar(clusters)]

Idents(seur_obj) <- "cluster_name"

# Split by cell type ------------------------------------------------------
for (cl in clusters) {
  if (!cl %in% unique(seur_obj$cluster_name)) {
    message("Skipping ", cl, " (not found in cluster_name)")
    next
  }

  cluster_obj <- subset(seur_obj, idents = cl)
  saveRDS(cluster_obj, file.path(data_root, paste0(cl, ".rds")))
}
