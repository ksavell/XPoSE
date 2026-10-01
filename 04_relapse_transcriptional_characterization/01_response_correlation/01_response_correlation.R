# ==============================================================================
# XPoSE-seq relapse response correlation
# Spearman correlation of RT Active vs Non-active transcriptional signatures
# within dmPFC, within vmPFC, and between regions.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# Configuration ----------------------------------------------------------
padj_threshold <- 0.05
lfc_threshold <- 0

analysis <- list(
  dmPFC = list(
    results_dir = "output/03_differential_expression/03_main_de/RT_active_RT_nonactive_dmPFC/results",
    clusters = c("ITL23", "ITL5", "ITL6", "CTL6", "ETL5", "NPL5", "Pvalb", "Sst")
  ),
  vmPFC = list(
    results_dir = "output/03_differential_expression/03_main_de/RT_active_RT_nonactive_vmPFC/results",
    clusters = c("ITL23", "ITL5", "ITL6", "ITvm", "CTL6", "ETL5", "NPL5", "Pvalb", "Sst")
  )
)

output_dir <- "output/04_relapse_transcriptional_characterization/01_response_correlation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Plot settings ----------------------------------------------------------
similarity_low <- "#2D8CB8"
similarity_mid <- "#FFFFFF"
similarity_high <- "#2D00FF"
diagonal_fill <- "#E5E5E5"
value_label_abs_threshold <- 0.01

# Read DE result tables --------------------------------------------------
read_signature <- function(region, cluster) {
  file <- file.path(
    analysis[[region]]$results_dir,
    paste0(cluster, "_DESeq2_results.csv")
  )

  if (!file.exists(file)) return(NULL)

  dat <- read_csv(file, show_col_types = FALSE)
  required <- c("gene", "log2FoldChange", "padj")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("Missing columns in ", file, ": ", paste(missing, collapse = ", "))
  }

  dat %>%
    transmute(
      gene = as.character(gene),
      log2FC = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    ) %>%
    filter(!duplicated(gene))
}

signatures <- list()
for (region in names(analysis)) {
  for (cluster in analysis[[region]]$clusters) {
    sig <- read_signature(region, cluster)
    if (!is.null(sig)) {
      signatures[[paste(region, cluster, sep = "__")]] <- sig
    }
  }
}

# Pairwise Spearman correlation -----------------------------------------
compare_pair <- function(region_a, cluster_a, region_b, cluster_b) {
  id_a <- paste(region_a, cluster_a, sep = "__")
  id_b <- paste(region_b, cluster_b, sep = "__")

  if (!id_a %in% names(signatures) || !id_b %in% names(signatures)) return(NULL)

  common <- inner_join(
    signatures[[id_a]] %>% rename(log2FC_a = log2FC, padj_a = padj),
    signatures[[id_b]] %>% rename(log2FC_b = log2FC, padj_b = padj),
    by = "gene"
  ) %>%
    filter(
      is.finite(log2FC_a),
      is.finite(log2FC_b),
      !is.na(padj_a),
      !is.na(padj_b)
    ) %>%
    mutate(
      sig_a = padj_a < padj_threshold & abs(log2FC_a) > lfc_threshold,
      sig_b = padj_b < padj_threshold & abs(log2FC_b) > lfc_threshold,
      in_union = sig_a | sig_b
    )

  union_dat <- common %>% filter(in_union)

  rho <- if (nrow(union_dat) >= 3) {
    suppressWarnings(cor(union_dat$log2FC_a, union_dat$log2FC_b, method = "spearman"))
  } else {
    NA_real_
  }

  p_value <- if (nrow(union_dat) >= 3) {
    suppressWarnings(
      cor.test(
        union_dat$log2FC_a,
        union_dat$log2FC_b,
        method = "spearman",
        exact = FALSE
      )$p.value
    )
  } else {
    NA_real_
  }

  tibble(
    region_a = region_a,
    cluster_a = cluster_a,
    region_b = region_b,
    cluster_b = cluster_b,
    n_tested_common = nrow(common),
    n_deg_a = sum(common$sig_a, na.rm = TRUE),
    n_deg_b = sum(common$sig_b, na.rm = TRUE),
    n_union = nrow(union_dat),
    spearman_rho = rho,
    spearman_p = p_value
  )
}

make_pair_table <- function(region_a, clusters_a, region_b, clusters_b) {
  rows <- list()
  k <- 1
  for (cluster_a in clusters_a) {
    for (cluster_b in clusters_b) {
      result <- compare_pair(region_a, cluster_a, region_b, cluster_b)
      if (!is.null(result)) {
        rows[[k]] <- result
        k <- k + 1
      }
    }
  }
  bind_rows(rows)
}

dm_pairs <- make_pair_table("dmPFC", analysis$dmPFC$clusters, "dmPFC", analysis$dmPFC$clusters)
vm_pairs <- make_pair_table("vmPFC", analysis$vmPFC$clusters, "vmPFC", analysis$vmPFC$clusters)
between_pairs <- make_pair_table("dmPFC", analysis$dmPFC$clusters, "vmPFC", analysis$vmPFC$clusters)

all_pairs <- bind_rows(
  dm_pairs %>% mutate(matrix = "dmPFC within", .before = 1),
  vm_pairs %>% mutate(matrix = "vmPFC within", .before = 1),
  between_pairs %>% mutate(matrix = "dmPFC x vmPFC", .before = 1)
)

write_csv(all_pairs, file.path(output_dir, "pairwise_spearman_statistics.csv"))
write_csv(dm_pairs, file.path(output_dir, "dmPFC_within_spearman.csv"))
write_csv(vm_pairs, file.path(output_dir, "vmPFC_within_spearman.csv"))
write_csv(between_pairs, file.path(output_dir, "dmPFC_x_vmPFC_spearman.csv"))

# Plotting ---------------------------------------------------------------
matrix_theme <- theme_minimal(base_family = "Arial", base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 7, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 8, face = "plain", colour = "black", hjust = 0),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(14, "mm"),
    legend.key.width = grid::unit(2.5, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  )

build_matrix_plot <- function(pair_data, row_clusters, col_clusters, title,
                              within_region = FALSE, outline_same_cluster = FALSE) {
  dat <- pair_data %>%
    mutate(
      row_cluster = factor(cluster_a, levels = rev(row_clusters)),
      col_cluster = factor(cluster_b, levels = col_clusters),
      diagonal = within_region & cluster_a == cluster_b,
      plot_value = if_else(diagonal, NA_real_, spearman_rho)
    )

  labels <- dat %>%
    filter(
      !diagonal,
      !is.na(spearman_rho),
      abs(spearman_rho) >= value_label_abs_threshold
    )

  p <- ggplot(dat, aes(x = col_cluster, y = row_cluster, fill = plot_value)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = similarity_low,
      mid = similarity_mid,
      high = similarity_high,
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      na.value = diagonal_fill,
      name = "Spearman\nrho"
    ) +
    geom_text(
      data = labels,
      aes(label = sprintf("%.2f", spearman_rho)),
      family = "Arial",
      size = 5.5 / ggplot2::.pt,
      colour = "black"
    ) +
    coord_fixed() +
    labs(title = title) +
    matrix_theme

  if (outline_same_cluster && !within_region) {
    p <- p +
      geom_tile(
        data = dat %>% filter(cluster_a == cluster_b),
        fill = NA,
        colour = "black",
        linewidth = 0.5
      )
  }

  p
}

p_dm <- build_matrix_plot(
  dm_pairs,
  analysis$dmPFC$clusters,
  analysis$dmPFC$clusters,
  "dmPFC",
  within_region = TRUE
)

p_vm <- build_matrix_plot(
  vm_pairs,
  analysis$vmPFC$clusters,
  analysis$vmPFC$clusters,
  "vmPFC",
  within_region = TRUE
)

p_between <- build_matrix_plot(
  between_pairs,
  analysis$dmPFC$clusters,
  analysis$vmPFC$clusters,
  "dmPFC x vmPFC",
  outline_same_cluster = TRUE
)

open_pdf <- function(file, width, height) {
  if (capabilities("aqua")) {
    quartz(type = "pdf", file = file, width = width, height = height,
           family = "Arial", pointsize = 7)
  } else {
    pdf(file = file, width = width, height = height,
        family = "Helvetica", pointsize = 7, useDingbats = FALSE)
  }
}

open_pdf(file.path(output_dir, "dmPFC_within_spearman.pdf"), 2.5, 2.35)
print(p_dm)
dev.off()

open_pdf(file.path(output_dir, "vmPFC_within_spearman.pdf"), 2.65, 2.5)
print(p_vm)
dev.off()

open_pdf(file.path(output_dir, "dmPFC_x_vmPFC_spearman.pdf"), 2.65, 2.5)
print(p_between)
dev.off()

combined <- (p_dm | p_vm | p_between) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

open_pdf(file.path(output_dir, "response_correlation_spearman_combined.pdf"), 7.1, 2.65)
print(combined)
dev.off()

message("Response-correlation analysis complete: ", output_dir)
