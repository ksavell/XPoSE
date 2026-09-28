# XPoSE-seq: cluster and annotate the main dataset
#
# Adaptation and consolidation for final repository: Katherine E. Savell
#
# Run from the base XPoSE repository directory.
# All four main experiences are clustered together across dmPFC and vmPFC.

rm(list = ls())

# Load packages -----------------------------------------------------------
library(Seurat)
library(ggplot2)

# Paths -------------------------------------------------------------------
input_file <- "output/01_metadata_clustering_qc/main_metadata_qc.rds"
config_file <- "01_metadata_clustering_qc/config/main_capture_config.csv"
function_file <- "01_metadata_clustering_qc/functions/clustering_functions.R"
output_dir <- "output/01_metadata_clustering_qc"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(function_file)

# Load MapMyCells annotations ---------------------------------------------
main_obj <- readRDS(input_file)
capture_config <- read_clustering_config(config_file)
main_obj <- join_mapmycells(main_obj, capture_config)

# Retain neuronal nuclei for clustering -----------------------------------
main_obj <- filter_for_neuronal_clustering(
  main_obj,
  min_features = 2000
)

# First clustering --------------------------------------------------------
main_obj <- cluster_first(
  main_obj,
  neigh_dim = 1:50,
  umap_dim = 1:50,
  resolution = 2
)

save_cluster_umap(
  main_obj,
  file.path(output_dir, "main_first_clustering_umap.pdf")
)

# Second clustering -------------------------------------------------------
# Hard-coded from manual review of the first clustering result.
main_clusters_keep <- as.character(0:41)

main_obj <- subset_recluster(
  main_obj,
  clusters_keep = main_clusters_keep,
  first_cluster_col = "RNA_snn_res.2",
  neigh_dim = 1:50,
  umap_dim = 1:40,
  resolution = 5
)

save_cluster_umap(
  main_obj,
  file.path(output_dir, "main_second_clustering_umap.pdf")
)

# Allen subclass annotation and manuscript cell-type labels ----------------
main_annotation <- annotate_clusters(main_obj)
main_obj <- main_annotation$object

write.csv(
  main_annotation$annotations,
  file = file.path(output_dir, "main_allen_annotations.csv"),
  row.names = FALSE
)

save_annotation_subclass_dotplot(
  main_obj,
  file.path(output_dir, "main_annotation_vs_subclass_dotplot.pdf")
)

saveRDS(
  main_obj,
  file = file.path(output_dir, "main_annotated.rds")
)
