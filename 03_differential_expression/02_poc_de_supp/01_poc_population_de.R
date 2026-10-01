#!/usr/bin/env Rscript
# XPoSE-seq: proof-of-concept population differential expression (Figure 5A)
#
# Run from the base XPoSE repository directory.
# Performs pseudobulk DESeq2 within each final cell population for the three
# population comparisons retained in Figure 5A.

suppressPackageStartupMessages({
  library(Seurat)
})

# Paths -------------------------------------------------------------------
input_file <- "output/01_metadata_clustering_qc/poc_combined_annotated.rds"
output_root <- "output/03_differential_expression/01_poc_population_de"

source("03_differential_expression/02_poc_de_sup/poc_de_functions.R")


dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Load POC object ---------------------------------------------------------
poc <- readRDS(input_file)

if (inherits(poc[["RNA"]], "Assay5")) {
  poc[["RNA"]] <- JoinLayers(poc[["RNA"]])
}

# Recreate the three analysis populations used for Figure 5A. Homecage nuclei
# are activity-agnostic; NC nuclei retain Active / Non-active assignments.
poc$de_population <- ifelse(
  poc$experience == "HC",
  "Homecage",
  ifelse(
    poc$population == "active",
    "Active",
    ifelse(poc$population == "non-active", "Non-active", NA_character_)
  )
)

poc <- subset(poc, subset = !is.na(de_population))
clusters <- sort(unique(poc$cluster_name))

# Figure 5A comparisons ---------------------------------------------------
comparisons <- list(
  Nonactive_vs_Homecage = c("de_population", "Non-active", "Homecage"),
  Active_vs_Homecage = c("de_population", "Active", "Homecage"),
  Active_vs_Nonactive = c("de_population", "Active", "Non-active")
)

for (comparison_name in names(comparisons)) {
  de_and_summary(
    seur_obj = poc,
    pair = comparisons[[comparison_name]],
    clusters = clusters,
    min_cell = 10,
    min_rat = 3,
    output_dir = file.path(output_root, comparison_name),
    save_dds = TRUE
  )
}
