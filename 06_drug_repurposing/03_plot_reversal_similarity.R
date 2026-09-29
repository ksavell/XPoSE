# ==============================================================================
# XPoSE-seq drug repurposing: selected-drug reversal similarity
# ==============================================================================
# Uses the RT ranked-signature gene-level export to make the two manuscript
# similarity matrices:
#   1. Binary reversed-gene Jaccard
#   2. Magnitude-weighted Jaccard using inverse rank strength
#
# The analytic unit is one cell type x human gene. Each drug pair is evaluated
# only over units represented for both drugs.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

cfg <- list(
  gene_level_file = "output/06_drug_repurposing/asgard/RT/selected_drug_ranked_signature_gene_level.csv.gz",
  output_dir = "output/06_drug_repurposing/figures",

  drug_order = c(
    "trazodone",
    "simvastatin",
    "noscapine",
    "dextromethorphan"
  ),

  short_labels = c(
    trazodone = "Traz.",
    simvastatin = "Sim.",
    noscapine = "Nosc.",
    dextromethorphan = "Dext."
  ),

  label_colors = c(
    trazodone = "#6FA6C4",
    simvastatin = "#78B642",
    noscapine = "#8B82D7",
    dextromethorphan = "#E68C96"
  ),

  heat_low = "#F2F2F2",
  heat_high = "#2D00FF",
  heat_limit = 0.75,
  diagonal_fill = "#E9E9E9",

  width = 2.15,
  height = 1.75
)

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cfg$gene_level_file)) {
  stop("Missing RT selected-drug gene-level file: ", cfg$gene_level_file)
}

gene_dat <- read_csv(cfg$gene_level_file, show_col_types = FALSE) %>%
  mutate(
    drug = tolower(trimws(as.character(drug))),
    cluster = as.character(cluster),
    human_gene = as.character(human_gene),
    reversed = as.logical(reversed),
    inverse_rank_strength = as.numeric(inverse_rank_strength),
    gene_key = paste(cluster, human_gene, sep = "|||")
  ) %>%
  filter(drug %in% cfg$drug_order) %>%
  group_by(drug, gene_key) %>%
  summarise(
    cluster = first(cluster),
    human_gene = first(human_gene),
    reversed = any(reversed %in% TRUE, na.rm = TRUE),
    inverse_rank_strength = max(inverse_rank_strength, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    inverse_rank_strength = if_else(
      reversed & is.finite(inverse_rank_strength),
      inverse_rank_strength,
      0
    )
  )

missing_drugs <- setdiff(cfg$drug_order, unique(gene_dat$drug))
if (length(missing_drugs) > 0) {
  stop("Selected drugs missing from RT gene-level export: ", paste(missing_drugs, collapse = ", "))
}

pairs <- combn(cfg$drug_order, 2, simplify = FALSE)
pairwise_rows <- list()

for (i in seq_along(pairs)) {
  drug_a <- pairs[[i]][1]
  drug_b <- pairs[[i]][2]

  a <- gene_dat %>%
    filter(drug == drug_a) %>%
    transmute(
      gene_key,
      reversed_a = reversed,
      weight_a = inverse_rank_strength
    )

  b <- gene_dat %>%
    filter(drug == drug_b) %>%
    transmute(
      gene_key,
      reversed_b = reversed,
      weight_b = inverse_rank_strength
    )

  common <- inner_join(a, b, by = "gene_key")

  rev_a <- common$reversed_a %in% TRUE
  rev_b <- common$reversed_b %in% TRUE

  n_shared <- sum(rev_a & rev_b)
  n_union <- sum(rev_a | rev_b)
  binary_jaccard <- if (n_union > 0) n_shared / n_union else NA_real_

  wa <- ifelse(rev_a & is.finite(common$weight_a), common$weight_a, 0)
  wb <- ifelse(rev_b & is.finite(common$weight_b), common$weight_b, 0)

  weighted_num <- sum(pmin(wa, wb), na.rm = TRUE)
  weighted_den <- sum(pmax(wa, wb), na.rm = TRUE)
  weighted_jaccard <- if (weighted_den > 0) weighted_num / weighted_den else NA_real_

  pairwise_rows[[i]] <- tibble(
    drug_a = drug_a,
    drug_b = drug_b,
    n_common_represented = nrow(common),
    n_reversed_a = sum(rev_a),
    n_reversed_b = sum(rev_b),
    n_shared_reversed = n_shared,
    n_reversed_union = n_union,
    binary_reversed_jaccard = binary_jaccard,
    weighted_jaccard_numerator = weighted_num,
    weighted_jaccard_denominator = weighted_den,
    inverse_rank_weighted_jaccard = weighted_jaccard
  )
}

pairwise <- bind_rows(pairwise_rows)
write_csv(pairwise, file.path(cfg$output_dir, "selected_drug_reversal_similarity_pairwise.csv"))

make_matrix <- function(value_col) {
  m <- matrix(
    NA_real_,
    nrow = length(cfg$drug_order),
    ncol = length(cfg$drug_order),
    dimnames = list(cfg$drug_order, cfg$drug_order)
  )

  for (i in seq_len(nrow(pairwise))) {
    a <- pairwise$drug_a[i]
    b <- pairwise$drug_b[i]
    value <- pairwise[[value_col]][i]
    m[a, b] <- value
    m[b, a] <- value
  }

  m
}

binary_matrix <- make_matrix("binary_reversed_jaccard")
weighted_matrix <- make_matrix("inverse_rank_weighted_jaccard")

write_csv(
  as.data.frame(binary_matrix, check.names = FALSE) %>% rownames_to_column("drug"),
  file.path(cfg$output_dir, "binary_reversal_jaccard_matrix.csv")
)
write_csv(
  as.data.frame(weighted_matrix, check.names = FALSE) %>% rownames_to_column("drug"),
  file.path(cfg$output_dir, "magnitude_weighted_jaccard_matrix.csv")
)

matrix_to_long <- function(m) {
  as.data.frame(as.table(m), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename(drug_y = Var1, drug_x = Var2, value = Freq) %>%
    mutate(
      drug_x = factor(drug_x, levels = cfg$drug_order),
      drug_y = factor(drug_y, levels = rev(cfg$drug_order)),
      diagonal = as.character(drug_x) == as.character(drug_y),
      value = if_else(diagonal, NA_real_, value),
      label = if_else(is.finite(value), as.character(round(100 * value)), "")
    )
}

pt_to_mm <- 25.4 / 72.27
outline_lwd <- 0.25 * pt_to_mm

make_heatmap <- function(m, title) {
  dat <- matrix_to_long(m)

  p <- ggplot(dat, aes(x = drug_x, y = drug_y)) +
    geom_tile(
      aes(fill = value),
      color = "white",
      linewidth = outline_lwd,
      width = 0.98,
      height = 0.98
    ) +
    geom_text(
      aes(label = label),
      family = "Arial",
      size = 7 / ggplot2::.pt,
      color = "black",
      na.rm = TRUE
    ) +
    scale_fill_gradient(
      low = cfg$heat_low,
      high = cfg$heat_high,
      limits = c(0, cfg$heat_limit),
      breaks = c(0, 0.25, 0.50, 0.75),
      labels = c("0", "25", "50", "75"),
      oob = scales::squish,
      na.value = cfg$diagonal_fill,
      name = "Jaccard index"
    ) +
    scale_x_discrete(labels = cfg$short_labels[cfg$drug_order], drop = FALSE) +
    scale_y_discrete(labels = cfg$short_labels[rev(cfg$drug_order)], drop = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_family = "Arial", base_size = 7) +
    theme(
      axis.line = element_line(linewidth = 0.5 * pt_to_mm, color = "black"),
      axis.ticks = element_blank(),
      axis.text.x = element_text(
        size = 7,
        family = "Arial",
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = unname(cfg$label_colors[cfg$drug_order])
      ),
      axis.text.y = element_text(
        size = 7,
        family = "Arial",
        color = unname(cfg$label_colors[rev(cfg$drug_order)])
      ),
      plot.title = element_text(size = 8, family = "Arial", hjust = 0, color = "black"),
      legend.title = element_text(size = 7, family = "Arial", color = "black", angle = 90),
      legend.text = element_text(size = 7, family = "Arial", color = "black"),
      legend.key.height = grid::unit(14, "mm"),
      legend.key.width = grid::unit(2.5, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )

  p
}

save_pdf <- function(plot_obj, filename) {
  out <- file.path(cfg$output_dir, filename)

  if (capabilities("aqua")) {
    quartz(
      file = out,
      type = "pdf",
      width = cfg$width,
      height = cfg$height,
      family = "Arial",
      pointsize = 7
    )
  } else {
    cairo_pdf(
      filename = out,
      width = cfg$width,
      height = cfg$height,
      family = "Arial",
      pointsize = 7
    )
  }

  print(plot_obj)
  dev.off()
  message("Saved: ", out)
}

save_pdf(
  make_heatmap(binary_matrix, "Reversed gene similarity"),
  "reversed_gene_similarity.pdf"
)

save_pdf(
  make_heatmap(weighted_matrix, "Reversal magnitude similarity"),
  "reversal_magnitude_similarity.pdf"
)
