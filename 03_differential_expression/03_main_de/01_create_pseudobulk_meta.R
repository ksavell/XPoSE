# XPoSE-seq: create pseudobulk count and metadata objects for main-dataset DE
#
# Run from the base XPoSE repository directory.
# Creates pseudobulk samples at the rat x cell-population x condition level for
# the combined dmPFC/vmPFC dataset, then saves region-specific subsets used by
# the configured DESeq2 contrasts.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

# Paths -------------------------------------------------------------------
input_file <- "output/01_metadata_clustering_qc/main_annotated.rds"
output_dir <- "output/03_differential_expression/03_main_de/pseudobulk"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Helper ------------------------------------------------------------------
build_pb_meta <- function(seur_obj, factors, label) {
  
  message("Building pseudobulk: ", label)
  
  factors_df <- FetchData(seur_obj, vars = factors) %>%
    mutate(across(everything(), as.character))
  
  # Aggregate on a temporary integer group rather than a concatenated text key.
  # This avoids ambiguity when factor levels contain spaces, hyphens, or underscores.
  combo_key <- do.call(
    paste,
    c(factors_df, list(sep = "\r"))
  )
  
  pb_group <- match(combo_key, unique(combo_key))
  seur_obj$pb_group <- pb_group
  
  # One metadata row per pseudobulk sample.
  meta <- factors_df %>%
    mutate(pb_group = pb_group) %>%
    distinct(pb_group, .keep_all = TRUE) %>%
    arrange(pb_group)
  
  # Human-readable sample identifier used to align metadata and count columns.
  sample_id <- do.call(
    paste,
    c(
      lapply(meta[factors], function(x) make.names(x)),
      list(sep = "_")
    )
  )
  
  if (anyDuplicated(sample_id)) {
    stop("Sanitized pseudobulk sample IDs are not unique for ", label, ".")
  }
  
  meta$sample_id <- sample_id
  meta$n_nuclei <- as.integer(table(factor(pb_group, levels = meta$pb_group)))
  
  # Sum raw counts across nuclei within each pseudobulk group.
  pb <- AggregateExpression(
    object = seur_obj,
    assays = "RNA",
    group.by = "pb_group",
    return.seurat = FALSE,
    normalization.method = "none",
    slot = "counts"
  )
  
  pb_counts <- pb$RNA
  
  # Seurat may prefix integer group names (for example, g1). Recover the
  # trailing integer and use it to map pseudobulk columns back to metadata.
  pb_group_from_counts <- suppressWarnings(
    as.integer(sub(".*?([0-9]+)$", "\\1", colnames(pb_counts)))
  )
  
  if (any(is.na(pb_group_from_counts))) {
    stop("Could not recover pseudobulk group IDs from count-matrix columns.")
  }
  
  common_groups <- pb_group_from_counts[
    pb_group_from_counts %in% meta$pb_group
  ]
  
  if (length(common_groups) == 0) {
    stop("No overlapping pseudobulk groups were found between counts and metadata.")
  }
  
  count_index <- match(common_groups, pb_group_from_counts)
  meta_index <- match(common_groups, meta$pb_group)
  
  pb_counts <- pb_counts[, count_index, drop = FALSE]
  meta <- meta[meta_index, , drop = FALSE]
  
  colnames(pb_counts) <- meta$sample_id
  rownames(meta) <- meta$sample_id
  
  meta <- meta %>%
    select(sample_id, all_of(factors), n_nuclei) %>%
    mutate(across(all_of(factors), factor)) %>%
    as.data.frame()
  
  rownames(meta) <- meta$sample_id
  
  if (!identical(colnames(pb_counts), rownames(meta))) {
    stop("Pseudobulk count columns and metadata rows are not aligned for ", label, ".")
  }
  
  list(
    pb = list(RNA = pb_counts),
    meta = meta
  )
}

# Load final main-dataset object -----------------------------------------
all <- readRDS(input_file)

factors <- c(
  "region",
  "population",
  "experience",
  "sex",
  "cluster_name",
  "ratID"
)

missing_factors <- setdiff(factors, colnames(all@meta.data))
if (length(missing_factors) > 0) {
  stop(
    "Required metadata columns are missing: ",
    paste(missing_factors, collapse = ", ")
  )
}

# Combined dmPFC + vmPFC pseudobulk -------------------------------------
# Region remains a grouping factor so paired dmPFC-vmPFC contrasts can be run.
res_split <- build_pb_meta(
  seur_obj = all,
  factors = factors,
  label = "main mPFC, region retained"
)

saveRDS(
  res_split$pb,
  file.path(output_dir, "main_pb_mPFC_split.rds")
)

saveRDS(
  res_split$meta,
  file.path(output_dir, "main_meta_mPFC_split.rds")
)

# Region-specific pseudobulk objects ------------------------------------
meta_dm <- res_split$meta %>%
  filter(region == "dmPFC")

pb_dm <- list(
  RNA = res_split$pb$RNA[, meta_dm$sample_id, drop = FALSE]
)

saveRDS(
  pb_dm,
  file.path(output_dir, "main_pb_dmPFC.rds")
)

saveRDS(
  meta_dm,
  file.path(output_dir, "main_meta_dmPFC.rds")
)

meta_vm <- res_split$meta %>%
  filter(region == "vmPFC")

pb_vm <- list(
  RNA = res_split$pb$RNA[, meta_vm$sample_id, drop = FALSE]
)

saveRDS(
  pb_vm,
  file.path(output_dir, "main_pb_vmPFC.rds")
)

saveRDS(
  meta_vm,
  file.path(output_dir, "main_meta_vmPFC.rds")
)
