# XPoSE-seq: assign metadata and remove doublets
#
# Original metadata workflow authors: Katherine E. Savell and Padmashri Saravanan
# Adaptation and consolidation for final repository: Katherine E. Savell
#
# Run from the base XPoSE repository directory.
# This script processes both the main and proof-of-concept datasets using
# dataset-specific capture configurations and metadata lookup tables.

rm(list = ls())

# Paths -------------------------------------------------------------------
function_file <- "01_metadata_clustering_qc/functions/metadata_functions.R"
output_dir <- "output/01_metadata_clustering_qc"

datasets <- data.frame(
  dataset = c("main", "poc"),
  capture_config_file = c(
    "01_metadata_clustering_qc/config/main_capture_config.csv",
    "01_metadata_clustering_qc/config/poc_capture_config.csv"
  ),
  metadata_lookup_file = c(
    "01_metadata_clustering_qc/config/main_metadata_lookup.csv",
    "01_metadata_clustering_qc/config/poc_metadata_lookup.csv"
  ),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Package checks -----------------------------------------------------------
required_packages <- c(
  "Seurat",
  "scDblFinder",
  "SingleCellExperiment"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

source(function_file)

# Process datasets ---------------------------------------------------------
for (dataset_i in seq_len(nrow(datasets))) {
  
  dataset_name <- datasets$dataset[dataset_i]
  message("Processing ", dataset_name, " dataset...")
  
  capture_config <- read_capture_config(
    datasets$capture_config_file[dataset_i]
  )
  metadata_lookup <- read_metadata_lookup(
    datasets$metadata_lookup_file[dataset_i]
  )
  
  config_captures <- capture_config$capture
  lookup_captures <- unique(metadata_lookup$capture)
  
  if (!setequal(config_captures, lookup_captures)) {
    stop(
      "Capture names differ between the ", dataset_name,
      " capture configuration and metadata lookup. Config-only: ",
      paste(setdiff(config_captures, lookup_captures), collapse = ", "),
      "; lookup-only: ",
      paste(setdiff(lookup_captures, config_captures), collapse = ", ")
    )
  }
  
  capture_objects <- setNames(
    vector("list", nrow(capture_config)),
    capture_config$capture
  )
  xpose_filter_summaries <- vector("list", nrow(capture_config))
  doublet_summaries <- vector("list", nrow(capture_config))
  
  # Process each nucleus capture ------------------------------------------
  for (i in seq_len(nrow(capture_config))) {
    
    capture_name <- capture_config$capture[i]
    region <- if ("region" %in% names(capture_config)) {
      capture_config$region[i]
    } else {
      NA_character_
    }
    
    message("  ", capture_name)
    
    valid_xpose_tags <- metadata_lookup$xpose_tag[
      metadata_lookup$capture == capture_name
    ]
    
    counts <- read_bd_table(capture_config$counts_file[i])
    tag_calls <- read_bd_table(capture_config$xpose_tag_calls_file[i])
    tag_reads <- read_bd_table(capture_config$xpose_tag_reads_file[i])
    
    seur_obj <- create_capture_seurat(
      counts = counts,
      tag_calls = tag_calls,
      tag_reads = tag_reads,
      capture = capture_name,
      region = region,
      valid_xpose_tags = valid_xpose_tags
    )
    
    seur_obj <- assign_experimental_metadata(
      seur_obj = seur_obj,
      metadata_lookup = metadata_lookup
    )
    
    # Remove between-XPoSE-tag multiplets, undetermined calls, and any tag
    # calls not represented in the lookup table for this capture.
    tag_filter <- filter_xpose_tag_calls(
      seur_obj = seur_obj,
      valid_xpose_tags = valid_xpose_tags
    )
    
    xpose_filter_summaries[[i]] <- tag_filter$summary
    seur_obj <- tag_filter$object
    
    # Detect and remove residual within-XPoSE-tag doublets.
    seur_obj <- run_within_xpose_doublet_detection(
      seur_obj = seur_obj,
      seed = 22
    )
    
    doublet_summaries[[i]] <- summarize_within_xpose_doublets(seur_obj)
    
    keep_cells <- colnames(seur_obj)[
      seur_obj$within_xpose_doublet == "singlet"
    ]
    seur_obj <- subset(seur_obj, cells = keep_cells)
    
    capture_objects[[capture_name]] <- seur_obj
    
    rm(counts, tag_calls, tag_reads, seur_obj, tag_filter)
    invisible(gc())
  }
  
  # Save filtering summaries ----------------------------------------------
  xpose_filter_summary <- do.call(rbind, xpose_filter_summaries)
  rownames(xpose_filter_summary) <- NULL
  
  write.csv(
    xpose_filter_summary,
    file = file.path(
      output_dir,
      paste0(dataset_name, "_xpose_tag_filter_summary.csv")
    ),
    row.names = FALSE
  )
  
  doublet_summary <- do.call(rbind, doublet_summaries)
  rownames(doublet_summary) <- NULL
  
  write.csv(
    doublet_summary,
    file = file.path(
      output_dir,
      paste0(dataset_name, "_within_xpose_doublet_summary.csv")
    ),
    row.names = FALSE
  )
  
  # Merge captures and save -----------------------------------------------
  dataset_seurat <- merge(
    x = capture_objects[[1]],
    y = capture_objects[-1],
    add.cell.ids = names(capture_objects),
    merge.data = FALSE,
    merge.dr = FALSE
  )
  
  saveRDS(
    dataset_seurat,
    file = file.path(
      output_dir,
      paste0(dataset_name, "_metadata_qc.rds")
    )
  )
  
  rm(
    capture_config,
    metadata_lookup,
    capture_objects,
    xpose_filter_summaries,
    doublet_summaries,
    xpose_filter_summary,
    doublet_summary,
    dataset_seurat
  )
  invisible(gc())
}