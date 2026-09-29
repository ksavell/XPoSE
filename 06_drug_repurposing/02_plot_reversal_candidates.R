# ==============================================================================
# XPoSE-seq drug repurposing: NC vs RT reversal-candidate comparison
# ==============================================================================
# Reproduces the manuscript scatter comparing final Asgard FDR values from the
# vmPFC NC Active-vs-Non-active and RT Active-vs-Non-active analyses.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

cfg <- list(
  asgard_root = "output/06_drug_repurposing/asgard",
  output_dir = "output/06_drug_repurposing/figures",
  fdr_threshold = 0.05,

  selected_drugs = c(
    "simvastatin",
    "trazodone",
    "noscapine",
    "dextromethorphan"
  ),

  selected_labels = c(
    simvastatin = "Simvastatin",
    trazodone = "Trazodone",
    noscapine = "Noscapine",
    dextromethorphan = "Dextromethorphan"
  ),

  category_colors = c(
    "ns" = "#D5D5D5",
    "NC selective" = "#3A8C87",
    "RT selective" = "#B64D79",
    "RT + NC" = "#6D637D"
  ),

  width = 3.15,
  height = 2.35,
  x_limit = c(0, 10),
  y_limit = c(0, 8)
)

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

read_scores <- function(experience) {
  file <- file.path(cfg$asgard_root, experience, "drug_scores.csv")
  if (!file.exists(file)) stop("Missing Asgard score table: ", file)

  read_csv(file, show_col_types = FALSE) %>%
    transmute(
      drug_key = tolower(trimws(as.character(drug))),
      therapeutic_score = as.numeric(therapeutic_score),
      p_value = as.numeric(p_value),
      fdr = as.numeric(fdr)
    )
}

rt <- read_scores("RT") %>%
  rename(
    therapeutic_score_RT = therapeutic_score,
    p_value_RT = p_value,
    fdr_RT = fdr
  )

nc <- read_scores("NC") %>%
  rename(
    therapeutic_score_NC = therapeutic_score,
    p_value_NC = p_value,
    fdr_NC = fdr
  )

plot_data <- full_join(rt, nc, by = "drug_key") %>%
  mutate(
    present_RT = !is.na(fdr_RT),
    present_NC = !is.na(fdr_NC),
    fdr_RT = replace_na(fdr_RT, 1),
    fdr_NC = replace_na(fdr_NC, 1),
    neglog10_fdr_RT = -log10(pmax(fdr_RT, .Machine$double.xmin)),
    neglog10_fdr_NC = -log10(pmax(fdr_NC, .Machine$double.xmin)),
    category = case_when(
      fdr_RT < cfg$fdr_threshold & fdr_NC < cfg$fdr_threshold ~ "RT + NC",
      fdr_RT < cfg$fdr_threshold ~ "RT selective",
      fdr_NC < cfg$fdr_threshold ~ "NC selective",
      TRUE ~ "ns"
    ),
    category = factor(
      category,
      levels = c("ns", "NC selective", "RT selective", "RT + NC")
    ),
    selected = drug_key %in% cfg$selected_drugs,
    label = if_else(
      selected,
      unname(cfg$selected_labels[drug_key]),
      NA_character_
    )
  )

write_csv(
  plot_data %>% arrange(category, desc(neglog10_fdr_RT), desc(neglog10_fdr_NC)),
  file.path(cfg$output_dir, "NC_vs_RT_reversal_candidates_plot_data.csv")
)

cutoff <- -log10(cfg$fdr_threshold)

pt_to_mm <- 25.4 / 72.27
axis_lwd <- 0.5 * pt_to_mm
outline_lwd <- 0.25 * pt_to_mm

p <- ggplot(plot_data, aes(x = neglog10_fdr_NC, y = neglog10_fdr_RT)) +
  geom_vline(
    xintercept = cutoff,
    linetype = "dotted",
    linewidth = axis_lwd,
    color = "black"
  ) +
  geom_hline(
    yintercept = cutoff,
    linetype = "dotted",
    linewidth = axis_lwd,
    color = "black"
  ) +
  geom_point(
    aes(fill = category),
    shape = 21,
    size = 1.8,
    stroke = outline_lwd,
    color = "white",
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = plot_data %>% filter(selected),
    aes(label = label),
    family = "Arial",
    size = 7 / ggplot2::.pt,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    segment.color = "#777777",
    segment.size = outline_lwd,
    seed = 4,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = cfg$category_colors,
    breaks = c("ns", "NC selective", "RT selective", "RT + NC"),
    drop = FALSE,
    name = NULL
  ) +
  scale_x_continuous(
    limits = cfg$x_limit,
    breaks = seq(cfg$x_limit[1], cfg$x_limit[2], by = 2),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    limits = cfg$y_limit,
    breaks = seq(cfg$y_limit[1], cfg$y_limit[2], by = 2),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "NC vs. RT reversal candidates",
    x = expression(-log[10](FDR)[NC]),
    y = expression(-log[10](FDR)[RT])
  ) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    axis.line = element_line(linewidth = axis_lwd, color = "black"),
    axis.ticks = element_line(linewidth = axis_lwd, color = "black"),
    axis.text = element_text(size = 7, family = "Arial", color = "black"),
    axis.title = element_text(size = 8, family = "Arial", color = "black"),
    plot.title = element_text(size = 8, family = "Arial", hjust = 0, color = "black"),
    legend.position = c(0.77, 0.76),
    legend.justification = c(0, 0.5),
    legend.text = element_text(size = 7, family = "Arial", color = "black"),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(3.2, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 2.1, color = NA)))

pdf_file <- file.path(cfg$output_dir, "NC_vs_RT_reversal_candidates.pdf")

if (capabilities("aqua")) {
  quartz(
    file = pdf_file,
    type = "pdf",
    width = cfg$width,
    height = cfg$height,
    family = "Arial",
    pointsize = 7
  )
} else {
  cairo_pdf(
    filename = pdf_file,
    width = cfg$width,
    height = cfg$height,
    family = "Arial",
    pointsize = 7
  )
}
print(p)
dev.off()

message("Saved: ", pdf_file)
