#!/usr/bin/env Rscript
# XPoSE-seq: plot POC active-fraction downsampling differential expression
#
# Generates the Figure 5 downsampling panel from de_summary_long.csv.
# Run locally on macOS for Quartz PDF output.
#
# Usage:
#   Rscript 05_plot_downsample.R --out_dir \
#     output/03_differential_expression/02_downsample_de/run_<run_id>

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option(
    "--out_dir",
    type = "character",
    help = "Collected downsampling output folder"
  )
)))

if (is.null(opt$out_dir)) {
  stop("--out_dir is required")
}

summary_file <- file.path(opt$out_dir, "de_summary_long.csv")
summary_df <- read_csv(summary_file, show_col_types = FALSE)

cluster_cols <- c(
  ITL23 = "#2EBF5E",
  ITL5  = "#50B2AD",
  ITL6  = "#58D2CF",
  CTL6  = "#2D8CB8",
  ETL5  = "#0D5A8B",
  Sst   = "#FF9900"
)

cluster_order <- c("ITL23", "ITL5", "ITL6", "CTL6", "ETL5", "Sst")
summary_df$cluster <- factor(summary_df$cluster, levels = cluster_order)

p <- ggplot(
  summary_df,
  aes(x = percentage, y = mean_DE, colour = cluster, group = cluster)
) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_DE - sd_DE),
      ymax = mean_DE + sd_DE
    ),
    width = 3,
    linewidth = 0.25,
    alpha = 0.75
  ) +
  geom_line(linewidth = 0.25, alpha = 0.4) +
  geom_point(size = 0.25, shape = 16, alpha = 0.75) +
  scale_colour_manual(
    values = cluster_cols,
    breaks = cluster_order,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 5, 25, 50, 75, 100),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dotted",
    colour = "black",
    linewidth = 0.25
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Active cells in NC pool (%)",
    y = "Mean DEGs",
    colour = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text = element_text(size = 7, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.key.width = grid::unit(10, "pt"),
    legend.key.height = grid::unit(7, "pt"),
    legend.spacing.x = grid::unit(2, "pt"),
    plot.title = element_blank()
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(linewidth = 0.5, size = 1.5)
    )
  )

quartz(
  type = "pdf",
  file = file.path(opt$out_dir, "de_by_active_fraction.pdf"),
  width = 2,
  height = 1.5,
  family = "Arial"
)
print(p)
dev.off()
