# XPoSE-seq: cluster and annotate the proof-of-concept dataset
#
# Adaptation and consolidation for final repository: Katherine E. Savell
#
# Run from the base XPoSE repository directory.
# The POC dataset is clustered as both an HC-only object and an HC+NC object.

rm(list = ls())

# Load packages -----------------------------------------------------------
library(Seurat)
library(ggplot2)

# Paths -------------------------------------------------------------------
input_file <- "output/01_metadata_clustering_qc/poc_metadata_qc.rds"
config_file <- "01_metadata_clustering_qc/config/poc_capture_config.csv"
function_file <- "01_metadata_clustering_qc/functions/clustering_functions.R"
output_dir <- "output/01_metadata_clustering_qc"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(function_file)

# Load MapMyCells annotations ---------------------------------------------
poc_obj <- readRDS(input_file)
capture_config <- read_clustering_config(config_file)
poc_obj <- join_mapmycells(poc_obj, capture_config)

# Neuron/non-neuron composition for the POC main figure ------------------
# Counts and percentages are calculated separately for each capture and rat
# before neuron or feature-count filtering.
poc_celltype_summary <- summarize_celltype_by_capture_rat(poc_obj)

write.csv(
  poc_celltype_summary,
  file = file.path(output_dir, "poc_neuron_non_neuron_by_capture_rat.csv"),
  row.names = FALSE
)

# Retain neuronal nuclei for clustering -----------------------------------
poc_obj <- filter_for_neuronal_clustering(
  poc_obj,
  min_features = 1000
)

# =========================================================================
# HC-only object
# =========================================================================
poc_hc <- subset(poc_obj, subset = experience == "HC")

# First clustering --------------------------------------------------------
poc_hc <- cluster_first(
  poc_hc,
  neigh_dim = 1:50,
  umap_dim = 1:50,
  resolution = 2
)

save_cluster_umap(
  poc_hc,
  file.path(output_dir, "poc_hc_first_clustering_umap.pdf")
)

# Second clustering -------------------------------------------------------
# Hard-coded from manual review; cluster 21 was Mbp+ and removed.
poc_hc_clusters_keep <- as.character(c(0:20, 22:23))

poc_hc <- subset_recluster(
  poc_hc,
  clusters_keep = poc_hc_clusters_keep,
  first_cluster_col = "RNA_snn_res.2",
  neigh_dim = 1:50,
  umap_dim = 1:40,
  resolution = 5
)

save_cluster_umap(
  poc_hc,
  file.path(output_dir, "poc_hc_second_clustering_umap.pdf")
)

# Allen subclass annotation -------------------------------------------------
poc_hc_annotation <- annotate_clusters(poc_hc)
poc_hc <- poc_hc_annotation$object

write.csv(
  poc_hc_annotation$annotations,
  file = file.path(output_dir, "poc_hc_allen_annotations.csv"),
  row.names = FALSE
)

save_annotation_subclass_dotplot(
  poc_hc,
  file.path(output_dir, "poc_hc_annotation_vs_subclass_dotplot.pdf")
)

saveRDS(
  poc_hc,
  file = file.path(output_dir, "poc_hc_annotated.rds")
)

# =========================================================================
# HC + NC combined object
# =========================================================================
poc_combined <- poc_obj

# First clustering --------------------------------------------------------
poc_combined <- cluster_first(
  poc_combined,
  neigh_dim = 1:50,
  umap_dim = 1:50,
  resolution = 2
)

save_cluster_umap(
  poc_combined,
  file.path(output_dir, "poc_combined_first_clustering_umap.pdf")
)

# Second clustering -------------------------------------------------------
# Hard-coded from manual review; cluster 24 was mural and 25 was Mbp+.
poc_combined_clusters_keep <- as.character(c(0:23, 26:27))

poc_combined <- subset_recluster(
  poc_combined,
  clusters_keep = poc_combined_clusters_keep,
  first_cluster_col = "RNA_snn_res.2",
  neigh_dim = 1:50,
  umap_dim = 1:40,
  resolution = 5
)

save_cluster_umap(
  poc_combined,
  file.path(output_dir, "poc_combined_second_clustering_umap.pdf")
)

# Allen subclass annotation -------------------------------------------------
poc_combined_annotation <- annotate_clusters(poc_combined)
poc_combined <- poc_combined_annotation$object

write.csv(
  poc_combined_annotation$annotations,
  file = file.path(output_dir, "poc_combined_allen_annotations.csv"),
  row.names = FALSE
)

save_annotation_subclass_dotplot(
  poc_combined,
  file.path(output_dir, "poc_combined_annotation_vs_subclass_dotplot.pdf")
)

saveRDS(
  poc_combined,
  file = file.path(output_dir, "poc_combined_annotated.rds")
)
