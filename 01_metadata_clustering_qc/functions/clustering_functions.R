# XPoSE-seq clustering and Allen MapMyCells annotation functions
#
# Consolidates the original cluster_first.R, subset_reclust.R, and
# MapMyCells annotation workflows.
# Adaptation and consolidation for final repository: Katherine E. Savell

# Cell classes retained for neuronal clustering.
neuron_classes_keep <- c(
  "01 IT-ET Glut",
  "02 NP-CT-L6b Glut",
  "06 CTX-CGE GABA",
  "07 CTX-MGE GABA",
  "08 CNU-MGE GABA"
)

# Final cell-type order used throughout the manuscript.
cluster_order <- c(
  "ITL23", "ITL5", "ITL6", "ITvm", "CTL6", "CTL6b", "ETL5", 
  "NPL5", "Pvalb", "Sst", "PvalbChand", "SstChodl", "Vip", "Lamp5", "Sncg"
)

# Allen subclass -> manuscript cell-type label.
subclass_collapse_map <- c(
  "007 L2/3 IT CTX Glut"       = "ITL23",
  "030 L6 CT CTX Glut"         = "CTL6",
  "010 IT AON-TT-DP Glut"      = "ITvm",
  "056 Sst Chodl Gaba"         = "SstChodl",
  "052 Pvalb Gaba"             = "Pvalb",
  "006 L4/5 IT CTX Glut"       = "ITL5",
  "004 L6 IT CTX Glut"         = "ITL6",
  "022 L5 ET CTX Glut"         = "ETL5",
  "053 Sst Gaba"               = "Sst",
  "029 L6b CTX Glut"           = "CTL6b",
  "032 L5 NP CTX Glut"         = "NPL5",
  "005 L5 IT CTX Glut"         = "ITL5",
  "046 Vip Gaba"               = "Vip",
  "049 Lamp5 Gaba"             = "Lamp5",
  "047 Sncg Gaba"              = "Sncg",
  "051 Pvalb chandelier Gaba"  = "PvalbChand",
  "050 Lamp5 Lhx6 Gaba"        = "Lamp5"
)

# Read capture-level MapMyCells result locations from the shared config.
read_clustering_config <- function(file) {
  config <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  
  required <- c("capture", "mapmycells_file")
  missing_cols <- setdiff(required, names(config))
  
  if (length(missing_cols) > 0) {
    stop(
      "Capture configuration is missing required clustering columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  config$capture <- trimws(config$capture)
  config$mapmycells_file <- trimws(config$mapmycells_file)
  config
}

# Join MapMyCells annotations to a Seurat object by the exported cell_id.
join_mapmycells <- function(seur_obj, capture_config) {
  annotation_cols <- c(
    "class_name",
    "subclass_name",
    "subclass_bootstrapping_probability",
    "supertype_name"
  )
  
  mapping_list <- lapply(seq_len(nrow(capture_config)), function(i) {
    mapping_file <- capture_config$mapmycells_file[i]
    capture_name <- capture_config$capture[i]
    
    if (!file.exists(mapping_file)) {
      stop("MapMyCells result not found: ", mapping_file)
    }
    
    mapping <- read.csv(
      mapping_file,
      comment.char = "#",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    required <- c("cell_id", annotation_cols)
    missing_cols <- setdiff(required, names(mapping))
    
    if (length(missing_cols) > 0) {
      stop(
        "MapMyCells file for ", capture_name,
        " is missing columns: ", paste(missing_cols, collapse = ", ")
      )
    }
    
    mapping$capture <- capture_name
    mapping[, c("cell_id", "capture", annotation_cols), drop = FALSE]
  })
  
  mapping <- do.call(rbind, mapping_list)
  
  if (anyDuplicated(mapping$cell_id)) {
    stop("Duplicate cell_id values were found across MapMyCells results.")
  }
  
  idx <- match(colnames(seur_obj), mapping$cell_id)
  
  for (annotation_col in annotation_cols) {
    seur_obj[[annotation_col]] <- mapping[[annotation_col]][idx]
  }
  
  class_number <- suppressWarnings(
    as.numeric(sub(" .*", "", seur_obj$class_name))
  )
  
  seur_obj$celltype <- ifelse(
    is.na(class_number),
    "unmapped",
    ifelse(class_number >= 30, "non-neuron", "neuron")
  )
  
  seur_obj
}

# Summarize POC neuron/non-neuron composition by capture and biological sample.
summarize_celltype_by_capture_rat <- function(seur_obj) {
  meta <- seur_obj@meta.data
  
  summary_df <- as.data.frame(
    table(
      capture = meta$capture,
      ratID = meta$ratID,
      celltype = meta$celltype,
      useNA = "ifany"
    ),
    stringsAsFactors = FALSE
  )
  
  names(summary_df)[names(summary_df) == "Freq"] <- "n_nuclei"
  summary_df <- summary_df[summary_df$n_nuclei > 0, , drop = FALSE]
  
  sample_total <- ave(
    summary_df$n_nuclei,
    summary_df$capture,
    summary_df$ratID,
    FUN = sum
  )
  
  summary_df$percent <- 100 * summary_df$n_nuclei / sample_total
  
  sample_key <- paste(meta$capture, meta$ratID, sep = "||")
  summary_key <- paste(summary_df$capture, summary_df$ratID, sep = "||")
  
  for (metadata_col in c("experience", "population", "sex")) {
    if (metadata_col %in% names(meta)) {
      lookup <- tapply(
        meta[[metadata_col]],
        sample_key,
        function(x) paste(unique(stats::na.omit(x)), collapse = ";")
      )
      summary_df[[metadata_col]] <- unname(lookup[summary_key])
    }
  }
  
  column_order <- c(
    "capture", "ratID", "experience", "population", "sex",
    "celltype", "n_nuclei", "percent"
  )
  column_order <- column_order[column_order %in% names(summary_df)]
  
  summary_df[, column_order, drop = FALSE]
}

# Retain neurons used for clustering and apply the dataset-specific feature cutoff.
filter_for_neuronal_clustering <- function(seur_obj,
                                           min_features,
                                           classes_keep = neuron_classes_keep) {
  seur_obj <- subset(seur_obj, subset = celltype == "neuron")
  seur_obj <- subset(seur_obj, subset = nFeature_RNA > min_features)
  seur_obj <- subset(seur_obj, subset = class_name %in% classes_keep)
  seur_obj
}

# Initial clustering, including the mitochondrial-content filter used originally.
cluster_first <- function(seur_obj,
                          neigh_dim = 1:50,
                          umap_dim = 1:50,
                          resolution = 2) {
  mitogenes <- c("ATP6", "COX1", "COX2", "CYTB", "ND1", "ND2", "ND4", "ND5", "ND6")
  seur_obj$percMito <- Seurat::PercentageFeatureSet(
    seur_obj,
    features = mitogenes
  )
  seur_obj <- subset(seur_obj, subset = percMito < 10)
  seur_obj$percMito <- NULL
  seur_obj <- Seurat::NormalizeData(seur_obj)
  seur_obj <- Seurat::FindVariableFeatures(
    seur_obj,
    selection.method = "vst",
    nfeatures = 2000
  )
  seur_obj <- Seurat::ScaleData(seur_obj, verbose = FALSE)
  seur_obj <- Seurat::RunPCA(
    seur_obj,
    npcs = max(c(neigh_dim, umap_dim)),
    verbose = FALSE
  )
  seur_obj <- Seurat::FindNeighbors(
    seur_obj,
    dims = neigh_dim,
    verbose = FALSE
  )
  seur_obj <- Seurat::RunUMAP(
    seur_obj,
    reduction = "pca",
    dims = umap_dim,
    verbose = FALSE
  )
  seur_obj <- Seurat::FindClusters(
    seur_obj,
    resolution = resolution,
    verbose = FALSE
  )
  
  seur_obj
}

# Remove manually reviewed first-pass clusters and recluster the retained nuclei.
subset_recluster <- function(seur_obj,
                             clusters_keep,
                             first_cluster_col = "RNA_snn_res.2",
                             neigh_dim = 1:50,
                             umap_dim = 1:40,
                             resolution = 5) {
  Seurat::Idents(seur_obj) <- first_cluster_col
  seur_obj <- subset(seur_obj, idents = clusters_keep)
  
  seur_obj <- Seurat::FindVariableFeatures(
    seur_obj,
    selection.method = "vst",
    nfeatures = 2000
  )
  seur_obj <- Seurat::ScaleData(seur_obj, verbose = FALSE)
  seur_obj <- Seurat::RunPCA(
    seur_obj,
    npcs = max(c(neigh_dim, umap_dim)),
    verbose = FALSE
  )
  seur_obj <- Seurat::FindNeighbors(
    seur_obj,
    dims = neigh_dim,
    verbose = FALSE
  )
  seur_obj <- Seurat::RunUMAP(
    seur_obj,
    reduction = "pca",
    dims = umap_dim,
    verbose = FALSE
  )
  seur_obj <- Seurat::FindClusters(
    seur_obj,
    resolution = resolution,
    verbose = FALSE
  )
  
  seur_obj
}

# Open a vector PDF using Quartz on macOS, with a standard PDF fallback.
open_vector_pdf <- function(file, width, height) {
  if (capabilities("aqua")) {
    grDevices::quartz(
      file = file,
      type = "pdf",
      width = width,
      height = height,
      family = "Arial",
      pointsize = 7
    )
  } else {
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      family = "Helvetica",
      pointsize = 7,
      useDingbats = FALSE
    )
  }
}

# Save a compact cluster-labeled UMAP for manual review.
save_cluster_umap <- function(seur_obj, file) {
  p <- Seurat::DimPlot(
    seur_obj,
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    Seurat::NoLegend() +
    ggplot2::theme_void(base_family = "Arial", base_size = 7)
  
  open_vector_pdf(file, width = 5, height = 5)
  print(p)
  grDevices::dev.off()
}

# Probability-weighted Allen label composition for each Seurat cluster.
weighted_annotation <- function(seur_obj,
                                cluster_col,
                                label_col,
                                probability_col) {
  meta <- seur_obj@meta.data[, c(cluster_col, label_col, probability_col), drop = FALSE]
  names(meta) <- c("cluster", "label", "probability")
  
  meta <- meta[
    !is.na(meta$cluster) &
      !is.na(meta$label) &
      !is.na(meta$probability),
    ,
    drop = FALSE
  ]
  
  weighted <- stats::aggregate(
    probability ~ cluster + label,
    data = meta,
    FUN = sum
  )
  
  clusters <- sort(unique(as.character(weighted$cluster)))
  
  annotation <- do.call(rbind, lapply(clusters, function(cluster_id) {
    x <- weighted[weighted$cluster == cluster_id, , drop = FALSE]
    x <- x[order(x$probability, decreasing = TRUE), , drop = FALSE]
    
    total_weight <- sum(x$probability)
    top1 <- x$probability[1] / total_weight
    top2 <- if (nrow(x) >= 2) x$probability[2] / total_weight else 0
    
    data.frame(
      cluster = cluster_id,
      auto_label = x$label[1],
      top_pct = round(top1 * 100, 1),
      margin_pct = round((top1 - top2) * 100, 1),
      near_tie = (top1 - top2) < 0.05,
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(annotation) <- NULL
  annotation
}

# Assign the dominant Allen subclass to each Seurat cluster and collapse it
# to the manuscript-facing cell-type label. Raw per-nucleus subclass_name and
# supertype_name fields from MapMyCells remain unchanged in the Seurat object.
annotate_clusters <- function(seur_obj,
                              cluster_col = "RNA_snn_res.5") {
  subclass_annotation <- weighted_annotation(
    seur_obj,
    cluster_col = cluster_col,
    label_col = "subclass_name",
    probability_col = "subclass_bootstrapping_probability"
  )
  
  cluster_ids <- as.character(seur_obj[[cluster_col, drop = TRUE]])
  
  subclass_lookup <- setNames(
    subclass_annotation$auto_label,
    subclass_annotation$cluster
  )
  
  allen_subclass <- unname(subclass_lookup[cluster_ids])
  
  collapsed_subclass <- unname(
    subclass_collapse_map[allen_subclass]
  )
  collapsed_subclass[is.na(collapsed_subclass)] <-
    allen_subclass[is.na(collapsed_subclass)]
  
  seur_obj$cluster_name <- collapsed_subclass
  Seurat::Idents(seur_obj) <- "cluster_name"
  
  annotation_table <- subclass_annotation
  names(annotation_table)[names(annotation_table) == "auto_label"] <-
    "allen_subclass_name"
  
  annotation_table$cluster_name <- unname(
    subclass_collapse_map[annotation_table$allen_subclass_name]
  )
  missing_label <- is.na(annotation_table$cluster_name)
  annotation_table$cluster_name[missing_label] <-
    annotation_table$allen_subclass_name[missing_label]
  
  list(
    object = seur_obj,
    annotations = annotation_table
  )
}

# Build the final annotation-vs-Allen-subclass dot plot.
make_annotation_subclass_dotplot <- function(seur_obj) {
  meta <- seur_obj@meta.data
  
  final_clusters <- cluster_order[cluster_order %in% unique(meta$cluster_name)]
  extra_clusters <- setdiff(unique(meta$cluster_name), final_clusters)
  final_clusters <- c(final_clusters, sort(extra_clusters))
  
  subclass_order <- names(subclass_collapse_map)[
    order(match(unname(subclass_collapse_map), cluster_order))
  ]
  subclass_order <- unique(subclass_order)
  subclass_order <- subclass_order[subclass_order %in% unique(meta$subclass_name)]
  subclass_order <- c(
    subclass_order,
    sort(setdiff(unique(meta$subclass_name), subclass_order))
  )
  
  plot_df <- do.call(rbind, lapply(final_clusters, function(cluster_name) {
    cluster_meta <- meta[meta$cluster_name == cluster_name, , drop = FALSE]
    
    do.call(rbind, lapply(subclass_order, function(subclass_name) {
      subclass_meta <- cluster_meta[
        cluster_meta$subclass_name == subclass_name,
        ,
        drop = FALSE
      ]
      
      data.frame(
        cluster_name = cluster_name,
        subclass_name = subclass_name,
        percent = 100 * nrow(subclass_meta) / nrow(cluster_meta),
        mean_probability = if (nrow(subclass_meta) > 0) {
          mean(
            subclass_meta$subclass_bootstrapping_probability,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  plot_df <- plot_df[plot_df$percent > 0, , drop = FALSE]
  plot_df$cluster_name <- factor(plot_df$cluster_name, levels = final_clusters)
  plot_df$subclass_name <- factor(plot_df$subclass_name, levels = rev(subclass_order))
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = cluster_name,
      y = subclass_name,
      size = percent,
      color = mean_probability
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(range = c(0.5, 4), name = "% nuclei") +
    ggplot2::scale_color_viridis_c(
      option = "magma",
      limits = c(0, 1),
      name = "Mean bootstrap\nprobability"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_classic(base_family = "Arial", base_size = 7) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.line = ggplot2::element_line(linewidth = 0.5),
      axis.ticks = ggplot2::element_line(linewidth = 0.5),
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 7)
    )
}

# Save the final annotation-vs-subclass dot plot.
save_annotation_subclass_dotplot <- function(seur_obj, file) {
  p <- make_annotation_subclass_dotplot(seur_obj)
  
  n_subclasses <- length(unique(seur_obj$subclass_name))
  plot_height <- max(3.5, 0.18 * n_subclasses + 1.5)
  
  open_vector_pdf(file, width = 6, height = plot_height)
  print(p)
  grDevices::dev.off()
}
