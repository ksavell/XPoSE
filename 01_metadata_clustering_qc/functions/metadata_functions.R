# XPoSE-seq metadata assignment and doublet-QC functions
#
# Original metadata workflow authors: Katherine E. Savell and Padmashri Saravanan
# Adaptation and consolidation for final repository: Katherine E. Savell

# Format numeric XPoSE-tag identifiers to the labels emitted by the BD pipeline.
format_xpose_tag <- function(x) {
  x <- trimws(as.character(x))

  numeric_tag <- grepl("^[0-9]+$", x)
  x[numeric_tag] <- sprintf("SampleTag%02d_mm", as.integer(x[numeric_tag]))

  short_bd_tag <- grepl("^SampleTag[0-9]_mm$", x)
  if (any(short_bd_tag)) {
    tag_number <- sub("^SampleTag([0-9])_mm$", "\\1", x[short_bd_tag])
    x[short_bd_tag] <- sprintf("SampleTag%02d_mm", as.integer(tag_number))
  }

  x
}

# Read and validate the capture configuration table.
read_capture_config <- function(file) {
  config <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "capture",
    "counts_file",
    "xpose_tag_calls_file",
    "xpose_tag_reads_file"
  )

  missing_cols <- setdiff(required, names(config))
  if (length(missing_cols) > 0) {
    stop(
      "Capture configuration is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  config$capture <- trimws(config$capture)
  if ("region" %in% names(config)) {
    config$region <- trimws(config$region)
  }

  if (anyDuplicated(config$capture)) {
    stop("Each capture must appear exactly once in the capture configuration.")
  }

  config
}

# Read and validate the experimental metadata lookup table.
read_metadata_lookup <- function(file) {
  lookup <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("capture", "xpose_tag")
  missing_cols <- setdiff(required, names(lookup))
  if (length(missing_cols) > 0) {
    stop(
      "Metadata lookup is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  lookup$capture <- trimws(lookup$capture)
  lookup$xpose_tag <- format_xpose_tag(lookup$xpose_tag)

  key <- paste(lookup$capture, lookup$xpose_tag, sep = "||")
  if (anyDuplicated(key)) {
    duplicate_keys <- unique(key[duplicated(key)])
    stop(
      "Metadata lookup contains duplicate capture/XPoSE-tag combinations: ",
      paste(duplicate_keys, collapse = ", ")
    )
  }

  lookup
}

# Read one BD Rhapsody CSV table.
read_bd_table <- function(file) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }

  read.csv(
    file,
    skip = 7,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = TRUE
  )
}

# Verify that the count, tag-call, and tag-read tables contain the same nuclei,
# then restore a common deterministic row order.
align_capture_tables <- function(counts, tag_calls, tag_reads, capture) {
  count_cells <- rownames(counts)
  call_cells <- rownames(tag_calls)
  read_cells <- rownames(tag_reads)

  if (!setequal(count_cells, call_cells) || !setequal(count_cells, read_cells)) {
    stop(
      "Cell IDs do not match across count/tag tables for capture ", capture, "."
    )
  }

  cell_order <- sort(count_cells)

  list(
    counts = counts[cell_order, , drop = FALSE],
    tag_calls = tag_calls[cell_order, , drop = FALSE],
    tag_reads = tag_reads[cell_order, , drop = FALSE]
  )
}

# Create one capture-level Seurat object and retain the raw XPoSE-tag calls.
create_capture_seurat <- function(counts,
                                  tag_calls,
                                  tag_reads,
                                  capture,
                                  region = NA_character_,
                                  valid_xpose_tags = character(0)) {
  aligned <- align_capture_tables(counts, tag_calls, tag_reads, capture)
  counts <- aligned$counts
  tag_calls <- aligned$tag_calls
  tag_reads <- aligned$tag_reads

  seur_obj <- Seurat::CreateSeuratObject(
    counts = t(as.data.frame(counts)),
    project = capture
  )

  seur_obj$raw_cell_id <- colnames(seur_obj)
  seur_obj$capture <- capture

  if (!is.na(region) && nzchar(region)) {
    seur_obj$region <- region
  }

  seur_obj$xpose_tag <- as.character(tag_calls[colnames(seur_obj), 1])
  seur_obj$percent_mt <- Seurat::PercentageFeatureSet(seur_obj, pattern = "^mt")

  # Add only the XPoSE-tag read-count columns that are valid for this capture.
  valid_xpose_tags <- unique(format_xpose_tag(valid_xpose_tags))

  for (xpose_tag in valid_xpose_tags) {
    tag_number <- sub("^SampleTag([0-9]{2})_mm$", "\\1", xpose_tag)
    raw_read_col <- paste0(xpose_tag, ".stAbO")
    metadata_col <- paste0("xpose_tag_", tag_number, "_reads")

    if (raw_read_col %in% colnames(tag_reads)) {
      seur_obj[[metadata_col]] <- tag_reads[colnames(seur_obj), raw_read_col]
    } else {
      warning(
        "XPoSE-tag read column ", raw_read_col,
        " was not found for capture ", capture, "."
      )
    }
  }

  seur_obj
}

# Assign all experimental metadata columns by matching capture + XPoSE-tag.
# Any call not represented in the lookup remains unassigned and is filtered
# before within-XPoSE-tag doublet detection.
assign_experimental_metadata <- function(seur_obj, metadata_lookup) {
  captures <- unique(seur_obj$capture)
  if (length(captures) != 1) {
    stop("assign_experimental_metadata() expects a single-capture Seurat object.")
  }

  capture <- captures[[1]]
  capture_lookup <- metadata_lookup[metadata_lookup$capture == capture, , drop = FALSE]

  if (nrow(capture_lookup) == 0) {
    stop("No metadata lookup entries were found for capture ", capture, ".")
  }

  idx <- match(seur_obj$xpose_tag, capture_lookup$xpose_tag)
  metadata_cols <- setdiff(names(capture_lookup), c("capture", "xpose_tag"))

  for (metadata_col in metadata_cols) {
    seur_obj[[metadata_col]] <- capture_lookup[[metadata_col]][idx]
  }

  seur_obj$metadata_assigned <- !is.na(idx)
  seur_obj
}

# Classify raw tag calls and remove between-XPoSE-tag multiplets,
# undetermined calls, and off-panel tag detections.
filter_xpose_tag_calls <- function(seur_obj, valid_xpose_tags) {
  valid_xpose_tags <- unique(format_xpose_tag(valid_xpose_tags))
  raw_calls <- as.character(seur_obj$xpose_tag)

  status <- rep("off_panel", length(raw_calls))
  status[raw_calls %in% valid_xpose_tags] <- "assigned"
  status[raw_calls == "Multiplet"] <- "multiplet"
  status[raw_calls == "Undetermined"] <- "undetermined"
  seur_obj$xpose_tag_status <- status

  summary_df <- as.data.frame(
    table(
      capture = seur_obj$capture,
      xpose_tag = seur_obj$xpose_tag,
      status = seur_obj$xpose_tag_status,
      useNA = "ifany"
    ),
    stringsAsFactors = FALSE
  )
  names(summary_df)[names(summary_df) == "Freq"] <- "n_nuclei"
  summary_df <- summary_df[summary_df$n_nuclei > 0, , drop = FALSE]

  keep_cells <- colnames(seur_obj)[seur_obj$xpose_tag_status == "assigned"]
  filtered_obj <- subset(seur_obj, cells = keep_cells)

  list(
    object = filtered_obj,
    summary = summary_df
  )
}

# Detect within-XPoSE-tag doublets separately within one capture.
run_within_xpose_doublet_detection <- function(seur_obj, seed = 42) {
  if (length(unique(seur_obj$capture)) != 1) {
    stop("Doublet detection must be run separately for each capture.")
  }

  if (any(seur_obj$xpose_tag_status != "assigned")) {
    stop(
      "Only assigned XPoSE-tag calls should enter within-tag doublet detection."
    )
  }

  sce <- Seurat::as.SingleCellExperiment(seur_obj)

  set.seed(seed)
  sce <- scDblFinder::scDblFinder(sce, samples = "xpose_tag")

  seur_obj$within_xpose_doublet <- sce$scDblFinder.class
  seur_obj$within_xpose_doublet_score <- sce$scDblFinder.score

  seur_obj
}

# Summarize within-XPoSE-tag doublet calls by capture and tag.
summarize_within_xpose_doublets <- function(seur_obj) {
  meta <- seur_obj@meta.data
  split_key <- interaction(meta$capture, meta$xpose_tag, drop = TRUE)
  split_meta <- split(meta, split_key)

  summaries <- lapply(split_meta, function(x) {
    n_nuclei <- nrow(x)
    n_doublets <- sum(x$within_xpose_doublet == "doublet", na.rm = TRUE)

    out <- data.frame(
      capture = unique(x$capture)[1],
      xpose_tag = unique(x$xpose_tag)[1],
      n_nuclei = n_nuclei,
      n_doublets = n_doublets,
      pct_doublets = round(100 * n_doublets / n_nuclei, 2),
      stringsAsFactors = FALSE
    )

    for (metadata_col in c("ratID", "orig_ratID", "experience", "sex", "group", "region")) {
      if (metadata_col %in% names(x)) {
        values <- unique(stats::na.omit(x[[metadata_col]]))
        out[[metadata_col]] <- if (length(values) == 1) values else paste(values, collapse = ";")
      }
    }

    out
  })

  do.call(rbind, summaries)
}

# Generate capture-level QC plots for scDblFinder calls.
make_doublet_qc_plot <- function(seur_obj, capture_label) {
  plot_obj <- seur_obj
  plot_obj <- Seurat::NormalizeData(plot_obj, verbose = FALSE)
  plot_obj <- Seurat::FindVariableFeatures(plot_obj, verbose = FALSE)
  plot_obj <- Seurat::ScaleData(plot_obj, verbose = FALSE)
  plot_obj <- Seurat::RunPCA(plot_obj, npcs = 20, verbose = FALSE)
  plot_obj <- Seurat::RunUMAP(plot_obj, dims = 1:20, verbose = FALSE)

  theme_qc <- ggplot2::theme_classic(base_family = "Arial", base_size = 7) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(size = 8, face = "plain"),
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 7),
      axis.line = ggplot2::element_line(linewidth = 0.5 / 2.845),
      axis.ticks = ggplot2::element_line(linewidth = 0.5 / 2.845)
    )

  score_plot <- Seurat::FeaturePlot(
    plot_obj,
    features = "within_xpose_doublet_score",
    pt.size = 0.25
  ) +
    ggplot2::scale_color_viridis_c(option = "magma") +
    ggplot2::labs(title = paste0(capture_label, " | doublet score")) +
    theme_qc

  call_plot <- Seurat::DimPlot(
    plot_obj,
    group.by = "within_xpose_doublet",
    pt.size = 0.25,
    cols = c(singlet = "grey80", doublet = "#D55E00")
  ) +
    ggplot2::labs(title = paste0(capture_label, " | doublet call")) +
    theme_qc

  meta <- seur_obj@meta.data
  bar_df <- as.data.frame(
    table(
      xpose_tag = meta$xpose_tag,
      call = meta$within_xpose_doublet
    ),
    stringsAsFactors = FALSE
  )

  bar_plot <- ggplot2::ggplot(
    bar_df,
    ggplot2::aes(x = xpose_tag, y = Freq, fill = call)
  ) +
    ggplot2::geom_col(position = "fill", linewidth = 0) +
    ggplot2::scale_fill_manual(
      values = c(singlet = "grey80", doublet = "#D55E00")
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      title = paste0(capture_label, " | doublet rate by XPoSE-tag"),
      x = "XPoSE-tag",
      y = "Proportion",
      fill = NULL
    ) +
    theme_qc +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  patchwork::wrap_plots(
    score_plot,
    call_plot,
    bar_plot,
    ncol = 2,
    widths = c(1, 1)
  )
}

# Open a vector PDF device. Quartz is used on macOS; standard PDF is used
# as a cross-platform fallback.
open_vector_pdf <- function(file, width, height) {
  if (identical(Sys.info()[["sysname"]], "Darwin") && capabilities("aqua")) {
    grDevices::quartz(
      type = "pdf",
      file = file,
      width = width,
      height = height
    )
  } else {
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      useDingbats = FALSE
    )
  }
}
