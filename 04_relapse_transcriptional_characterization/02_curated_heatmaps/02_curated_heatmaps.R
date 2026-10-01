# ==============================================================================
# XPoSE-seq relapse curated DEG heatmaps
# Joins final curated gene/category assignments directly to current DE tables.
# ==============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readr)
  library(grid)
})

# Configuration ----------------------------------------------------------
reference_file <- "04_relapse_transcriptional_characterization/02_curated_heatmaps/curated_heatmap_reference.csv"

de_results_dirs <- c(
  dmPFC = "output/03_differential_expression/03_main_de/RT_active_RT_nonactive_dmPFC/results",
  vmPFC = "output/03_differential_expression/03_main_de/RT_active_RT_nonactive_vmPFC/results"
)

output_dir <- "output/04_relapse_transcriptional_characterization/02_curated_heatmaps/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Choose one: main_compact, supplemental_up_only, supplemental_all
panel_set <- "main_compact"

write_region_heatmaps <- TRUE
write_paired_heatmap <- TRUE
paired_gene_set_mode <- "union"  # union or intersection

padj_threshold <- 0.05

main_category_order <- c(
  "Activity-dependent IEG",
  "Glutamatergic",
  "GABA",
  "Neuromodulatory receptors",
  "Calcium channels",
  "Potassium channels",
  "Sodium channels",
  "Neuropeptides",
  "Synaptic plasticity",
  "Epigenetic regulators"
)

supplemental_extra_categories <- c(
  "Sigma receptor system",
  "Related GPCRs",
  "Chloride channels"
)

category_order <- if (panel_set == "main_compact") {
  main_category_order
} else {
  c(main_category_order, supplemental_extra_categories)
}

category_display <- c(
  "Activity-dependent IEG" = "Activity\ndependent",
  "Glutamatergic" = "Glutamatergic",
  "GABA" = "GABA",
  "Neuromodulatory receptors" = "Neuromodulatory\nreceptors",
  "Calcium channels" = "Calcium channels",
  "Potassium channels" = "Potassium channels",
  "Sodium channels" = "Sodium channels",
  "Neuropeptides" = "Neuropeptides",
  "Synaptic plasticity" = "Synaptic\nplasticity",
  "Epigenetic regulators" = "Epigenetic\nregulators",
  "Sigma receptor system" = "Sigma receptor\nsystem",
  "Related GPCRs" = "Related GPCRs",
  "Chloride channels" = "Chloride channels"
)

cluster_order <- list(
  dmPFC = c("ITL23", "ITL5", "ITL6", "CTL6", "ETL5", "NPL5", "Pvalb", "Sst"),
  vmPFC = c("ITL23", "ITL5", "ITL6", "ITvm", "CTL6", "ETL5", "NPL5", "Pvalb", "Sst")
)

cluster_display <- c(
  ITL23 = "IT L2/3", ITL5 = "IT L5", ITL6 = "IT L6", ITvm = "IT vm",
  CTL6 = "CT L6", CTL6b = "CT L6b", ETL5 = "ET L5", NPL5 = "NP L5",
  Pvalb = "Pvalb", PvalbChand = "PvalbChand", Sst = "Sst",
  SstChodl = "SstChodl", Sncg = "Sncg", Vip = "Vip", Lamp5 = "Lamp5"
)

cluster_colors <- c(
  ITL23 = "#2EBF5E", ITL5 = "#50B2AD", ITL6 = "#58D2CF", ITvm = "#B1DE7D",
  CTL6 = "#2D8CB8", CTL6b = "#7044AA", ETL5 = "#0D5A8B", NPL5 = "#3E9E64",
  Pvalb = "#B9342C", PvalbChand = "#FF2D4E", Sst = "#FF9900",
  SstChodl = "#B1B10C", Sncg = "#D3408D", Vip = "#B864CC", Lamp5 = "#DA808C"
)

lfc_min <- -2.5
lfc_mid <- 0
lfc_max <- 2.5
lfc_colors <- c("#444444", "#FFFFFF", "#000080")
significance_box_col <- "#8B174D"
significance_box_lwd <- 0.55

region_width_in <- 3.25
region_height_in <- 5.0
paired_width_in <- 6.0
paired_height_in <- 5.0

# Load curated reference -------------------------------------------------
reference <- read_csv(reference_file, show_col_types = FALSE)

panel_specs <- list(
  main_compact = c(include = "main_compact", rank = "main_display_rank"),
  supplemental_up_only = c(include = "supplemental_up_only", rank = "supplemental_up_display_rank"),
  supplemental_all = c(include = "supplemental_all", rank = "supplemental_all_display_rank")
)

if (!panel_set %in% names(panel_specs)) {
  stop("panel_set must be one of: ", paste(names(panel_specs), collapse = ", "))
}

include_col <- panel_specs[[panel_set]][["include"]]
rank_col <- panel_specs[[panel_set]][["rank"]]

reference_use <- reference %>%
  filter(.data[[include_col]], final_category %in% category_order) %>%
  mutate(
    display_rank = as.numeric(.data[[rank_col]]),
    final_category = factor(final_category, levels = category_order)
  ) %>%
  arrange(region, final_category, display_rank, gene)

# Read DE tables ---------------------------------------------------------
read_region_de <- function(region) {
  results_dir <- de_results_dirs[[region]]
  files <- list.files(results_dir, pattern = "_DESeq2_results\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No DE result files found in: ", results_dir)

  bind_rows(lapply(files, function(file) {
    cluster <- sub("_DESeq2_results\\.csv$", "", basename(file))
    read_csv(file, show_col_types = FALSE) %>%
      transmute(
        region = region,
        cluster = cluster,
        gene = as.character(gene),
        log2FC = as.numeric(log2FoldChange),
        padj = as.numeric(padj),
        significant = !is.na(padj) & padj < padj_threshold
      )
  }))
}

de_long <- bind_rows(
  read_region_de("dmPFC"),
  read_region_de("vmPFC")
)

# Build heatmap matrices -------------------------------------------------
build_region_data <- function(region, row_meta = NULL) {
  if (is.null(row_meta)) {
    row_meta <- reference_use %>%
      filter(region == !!region) %>%
      transmute(
        final_category = as.character(final_category),
        gene,
        display_rank
      )
  }

  clusters <- cluster_order[[region]]
  dat <- de_long %>%
    filter(region == !!region, cluster %in% clusters)

  mat <- matrix(
    NA_real_,
    nrow = nrow(row_meta),
    ncol = length(clusters),
    dimnames = list(row_meta$gene, clusters)
  )

  sig <- matrix(
    FALSE,
    nrow = nrow(row_meta),
    ncol = length(clusters),
    dimnames = list(row_meta$gene, clusters)
  )

  for (i in seq_len(nrow(row_meta))) {
    for (j in seq_along(clusters)) {
      hit <- dat %>%
        filter(gene == row_meta$gene[i], cluster == clusters[j]) %>%
        slice_head(n = 1)

      if (nrow(hit) == 1) {
        mat[i, j] <- hit$log2FC
        sig[i, j] <- hit$significant
      }
    }
  }

  list(region = region, mat = mat, sig = sig, row_meta = row_meta, clusters = clusters)
}

make_heatmap <- function(x, title = NULL, show_row_names = TRUE) {
  row_split <- factor(
    x$row_meta$final_category,
    levels = category_order,
    labels = unname(category_display[category_order])
  )

  display_names <- unname(cluster_display[x$clusters])
  display_names[is.na(display_names)] <- x$clusters[is.na(display_names)]
  colnames(x$mat) <- display_names
  colnames(x$sig) <- display_names

  column_text_colors <- unname(cluster_colors[x$clusters])
  column_text_colors[is.na(column_text_colors)] <- "#000000"

  col_fun <- circlize::colorRamp2(c(lfc_min, lfc_mid, lfc_max), lfc_colors)

  Heatmap(
    x$mat,
    name = "log2FC",
    col = col_fun,
    na_col = "#F2F2F2",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = row_split,
    cluster_row_slices = FALSE,
    row_gap = unit(1, "mm"),
    column_gap = unit(0.25, "mm"),
    rect_gp = gpar(col = NA),
    show_row_names = show_row_names,
    row_names_side = "left",
    row_title_side = "left",
    row_names_gp = gpar(fontsize = 5.5, fontfamily = "Arial", fontface = "italic"),
    row_title_gp = gpar(fontsize = 7, fontfamily = "Arial"),
    row_title_rot = 0,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_centered = FALSE,
    column_names_gp = gpar(fontsize = 7, fontfamily = "Arial", col = column_text_colors),
    column_title = title,
    column_title_gp = gpar(fontsize = 7, fontfamily = "Arial"),
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x_pos, y_pos, width, height, fill) {
      if (isTRUE(x$sig[i, j])) {
        grid.rect(
          x = x_pos, y = y_pos, width = width, height = height,
          gp = gpar(fill = NA, col = significance_box_col, lwd = significance_box_lwd)
        )
      }
    }
  )
}

open_pdf <- function(file, width, height) {
  if (capabilities("aqua")) {
    quartz(type = "pdf", file = file, width = width, height = height,
           family = "Arial", pointsize = 7)
  } else {
    pdf(file = file, width = width, height = height,
        family = "Helvetica", pointsize = 7, useDingbats = FALSE)
  }
}

write_plot_data <- function(prefix, x) {
  long <- expand.grid(
    gene = rownames(x$mat),
    cluster = x$clusters,
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      x$row_meta %>% select(gene, final_category, display_rank),
      by = "gene"
    )

  long$log2FC <- as.vector(x$mat)
  long$significant <- as.vector(x$sig)
  long$region <- x$region

  write_csv(
    long %>% select(region, final_category, display_rank, gene, cluster, log2FC, significant),
    file.path(output_dir, paste0(prefix, "_heatmap_values.csv"))
  )
}

# Standalone log2FC legend ----------------------------------------------
legend_file <- file.path(output_dir, "log2FC_color_legend.pdf")
col_fun <- circlize::colorRamp2(c(lfc_min, lfc_mid, lfc_max), lfc_colors)
lgd <- Legend(
  title = "log2FC",
  col_fun = col_fun,
  at = c(lfc_min, 0, lfc_max),
  labels = c(as.character(lfc_min), "0", as.character(lfc_max)),
  direction = "vertical",
  legend_height = unit(0.85, "in"),
  grid_width = unit(0.12, "in"),
  border = "#666666",
  title_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  labels_gp = gpar(fontsize = 6.5, fontfamily = "Arial")
)
open_pdf(legend_file, 0.85, 1.65)
grid.newpage()
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center", "center"))
dev.off()

# Region heatmaps --------------------------------------------------------
region_data <- list()

if (write_region_heatmaps) {
  for (region in c("dmPFC", "vmPFC")) {
    x <- build_region_data(region)
    region_data[[region]] <- x

    ht <- make_heatmap(x)
    file <- file.path(output_dir, paste0(region, "_", panel_set, "_heatmap.pdf"))
    open_pdf(file, region_width_in, region_height_in)
    draw(ht, padding = unit(c(1.5, 1.5, 1.5, 1.5), "mm"))
    dev.off()

    write_plot_data(paste0(region, "_", panel_set), x)
  }
}

# Paired-region heatmap --------------------------------------------------
if (write_paired_heatmap) {
  dm_ref <- reference_use %>% filter(region == "dmPFC")
  vm_ref <- reference_use %>% filter(region == "vmPFC")

  dm_keys <- paste(dm_ref$final_category, dm_ref$gene, sep = "||")
  vm_keys <- paste(vm_ref$final_category, vm_ref$gene, sep = "||")

  keep_keys <- if (paired_gene_set_mode == "intersection") {
    intersect(dm_keys, vm_keys)
  } else {
    union(dm_keys, vm_keys)
  }

  combined_ref <- bind_rows(dm_ref, vm_ref) %>%
    mutate(row_key = paste(final_category, gene, sep = "||")) %>%
    filter(row_key %in% keep_keys) %>%
    group_by(row_key, final_category, gene) %>%
    summarise(display_rank = min(display_rank, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      final_category = factor(final_category, levels = category_order),
      category_rank = match(final_category, category_order)
    ) %>%
    arrange(category_rank, display_rank, gene) %>%
    transmute(final_category = as.character(final_category), gene, display_rank)

  dm_pair <- build_region_data("dmPFC", row_meta = combined_ref)
  vm_pair <- build_region_data("vmPFC", row_meta = combined_ref)

  ht_dm <- make_heatmap(dm_pair, title = "dmPFC", show_row_names = TRUE)
  ht_vm <- make_heatmap(vm_pair, title = "vmPFC", show_row_names = FALSE)

  file <- file.path(
    output_dir,
    paste0("paired_dmPFC_vmPFC_", panel_set, "_", paired_gene_set_mode, "_heatmap.pdf")
  )

  open_pdf(file, paired_width_in, paired_height_in)
  draw(ht_dm + ht_vm, padding = unit(c(1.5, 1.5, 1.5, 1.5), "mm"))
  dev.off()

  write_plot_data(paste0("paired_", panel_set, "_dmPFC"), dm_pair)
  write_plot_data(paste0("paired_", panel_set, "_vmPFC"), vm_pair)
}

message("Curated heatmaps complete: ", output_dir)
