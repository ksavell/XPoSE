# XPoSE-seq proof-of-concept dataset: assign metadata and remove doublets
#
# Original metadata workflow authors: Katherine E. Savell and Padmashri Saravanan
# Adaptation and consolidation for final repository: Katherine E. Savell
#
# Run from the base XPoSE repository directory.
# The proof-of-concept dataset contains two nucleus captures and intentionally
# uses XPoSE-tags 2-9 only; tag calls outside that set are removed because they
# are not represented in the metadata lookup table.

rm(list = ls())

# Paths -------------------------------------------------------------------
function_file <- "01_metadata_clustering_qc/functions/metadata_functions.R"
capture_config_file <- "01_metadata_clustering_qc/config/poc_capture_config.csv"
metadata_lookup_file <- "01_metadata_clustering_qc/config/poc_metadata_lookup.csv"
output_dir <- "output/01_metadata_clustering_qc"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Package checks -----------------------------------------------------------
required_packages <- c(
  "Seurat",
  "scDblFinder",
  "SingleCellExperiment",
  "ggplot2",
  "patchwork",
  "scales"
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

# Load configuration and metadata lookup ----------------------------------
capture_config <- read_capture_config(capture_config_file)
metadata_lookup <- read_metadata_lookup(metadata_lookup_file)

config_captures <- capture_config$capture
lookup_captures <- unique(metadata_lookup$capture)

if (!setequal(config_captures, lookup_captures)) {
  stop(
    "Capture names differ between the capture configuration and metadata lookup. ",
    "Config-only: ", paste(setdiff(config_captures, lookup_captures), collapse = ", "),
    "; lookup-only: ", paste(setdiff(lookup_captures, config_captures), collapse = ", ")
  )
}

# Process each capture -----------------------------------------------------
capture_objects <- setNames(vector("list", nrow(capture_config)), capture_config$capture)
xpose_filter_summaries <- vector("list", nrow(capture_config))
doublet_summaries <- vector("list", nrow(capture_config))
doublet_qc_plots <- vector("list", nrow(capture_config))

for (i in seq_len(nrow(capture_config))) {
  capture <- capture_config$capture[i]
  region <- if ("region" %in% names(capture_config)) capture_config$region[i] else NA_character_

  message("Processing ", capture, "...")

  valid_xpose_tags <- metadata_lookup$xpose_tag[
    metadata_lookup$capture == capture
  ]

  counts <- read_bd_table(capture_config$counts_file[i])
  tag_calls <- read_bd_table(capture_config$xpose_tag_calls_file[i])
  tag_reads <- read_bd_table(capture_config$xpose_tag_reads_file[i])

  seur_obj <- create_capture_seurat(
    counts = counts,
    tag_calls = tag_calls,
    tag_reads = tag_reads,
    capture = capture,
    region = region,
    valid_xpose_tags = valid_xpose_tags
  )

  seur_obj <- assign_experimental_metadata(
    seur_obj = seur_obj,
    metadata_lookup = metadata_lookup
  )

  # This removes Multiplet/Undetermined calls as well as the occasional
  # off-panel SampleTag10_mm/SampleTag12_mm detections from the POC dataset.
  tag_filter <- filter_xpose_tag_calls(
    seur_obj = seur_obj,
    valid_xpose_tags = valid_xpose_tags
  )

  xpose_filter_summaries[[i]] <- tag_filter$summary
  seur_obj <- tag_filter$object

  message(
    capture, ": ", ncol(seur_obj),
    " nuclei remain after XPoSE-tag assignment filtering."
  )

  # Detect residual within-XPoSE-tag doublets. No singleton-tag workaround is
  # used; only valid XPoSE-tag assignments enter scDblFinder.
  seur_obj <- run_within_xpose_doublet_detection(
    seur_obj = seur_obj,
    seed = 42
  )

  doublet_summaries[[i]] <- summarize_within_xpose_doublets(seur_obj)
  doublet_qc_plots[[i]] <- make_doublet_qc_plot(seur_obj, capture)

  keep_cells <- colnames(seur_obj)[
    seur_obj$within_xpose_doublet == "singlet"
  ]
  seur_obj <- subset(seur_obj, cells = keep_cells)

  message(
    capture, ": ", ncol(seur_obj),
    " nuclei remain after within-XPoSE-tag doublet removal."
  )

  capture_objects[[capture]] <- seur_obj

  rm(counts, tag_calls, tag_reads, seur_obj)
  invisible(gc())
}

# Save filtering summaries -------------------------------------------------
xpose_filter_summary <- do.call(rbind, xpose_filter_summaries)
rownames(xpose_filter_summary) <- NULL

write.csv(
  xpose_filter_summary,
  file = file.path(output_dir, "poc_xpose_tag_filter_summary.csv"),
  row.names = FALSE
)

doublet_summary <- do.call(rbind, doublet_summaries)
rownames(doublet_summary) <- NULL

write.csv(
  doublet_summary,
  file = file.path(output_dir, "poc_within_xpose_doublet_summary.csv"),
  row.names = FALSE
)

# Save doublet QC plots ----------------------------------------------------
qc_pdf <- file.path(output_dir, "poc_within_xpose_doublet_qc.pdf")
open_vector_pdf(qc_pdf, width = 7.2, height = 5.2)

for (i in seq_along(doublet_qc_plots)) {
  print(doublet_qc_plots[[i]])
}

cross_capture_plot <- ggplot2::ggplot(
  doublet_summary,
  ggplot2::aes(x = xpose_tag, y = pct_doublets / 100)
) +
  ggplot2::geom_col(fill = "grey55", linewidth = 0) +
  ggplot2::facet_wrap(~capture, ncol = 2) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::labs(
    x = "XPoSE-tag",
    y = "Within-tag doublet rate"
  ) +
  ggplot2::theme_classic(base_family = "Arial", base_size = 7) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 8),
    axis.text = ggplot2::element_text(size = 7),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 7),
    axis.line = ggplot2::element_line(linewidth = 0.5 / 2.845),
    axis.ticks = ggplot2::element_line(linewidth = 0.5 / 2.845)
  )

print(cross_capture_plot)
grDevices::dev.off()

# Merge both captures ------------------------------------------------------
# Capture prefixes make cell names unique after merging; raw_cell_id retains
# the original BD nucleus barcode for traceability.
poc_seurat <- merge(
  x = capture_objects[[1]],
  y = capture_objects[-1],
  add.cell.ids = names(capture_objects),
  merge.data = FALSE,
  merge.dr = FALSE
)

saveRDS(
  poc_seurat,
  file = file.path(output_dir, "poc_metadata_qc.rds")
)
