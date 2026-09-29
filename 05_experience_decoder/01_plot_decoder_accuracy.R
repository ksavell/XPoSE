# ==============================================================================
# XPoSE-seq experience decoder: regional percent correct
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

input_file <- "output/05_experience_decoder/decoder_percent_correct.csv"
output_file <- "output/05_experience_decoder/decoder_percent_correct.pdf"

region_order <- c("dmPFC", "vmPFC")
region_colors <- c(
  dmPFC = "#8B3FA8",
  vmPFC = "#F07800"
)

plot_width_in <- 1.35
plot_height_in <- 1.50

pt <- function(x) x / 2.835

plot_data <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(region %in% region_order) %>%
  mutate(region = factor(region, levels = region_order)) %>%
  arrange(region)

if (nrow(plot_data) != 2L || anyDuplicated(plot_data$region)) {
  stop("Expected exactly one dmPFC row and one vmPFC row in ", input_file)
}

p <- ggplot(
  plot_data,
  aes(x = region, y = percent_correct, fill = region)
) +
  geom_hline(
    yintercept = 50,
    linewidth = pt(0.25),
    linetype = "dashed",
    colour = "grey60"
  ) +
  geom_col(
    width = 0.56,
    colour = "black",
    linewidth = pt(0.25)
  ) +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "% correct"
  ) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    axis.text = element_text(
      family = "Arial",
      size = 7,
      colour = "black"
    ),
    axis.title = element_text(
      family = "Arial",
      size = 8,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = pt(0.5),
      colour = "black"
    ),
    axis.ticks = element_line(
      linewidth = pt(0.5),
      colour = "black"
    ),
    axis.ticks.length = grid::unit(1.5, "pt"),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  )

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

if (capabilities("aqua")) {
  quartz(
    type = "pdf",
    file = output_file,
    width = plot_width_in,
    height = plot_height_in,
    family = "Arial"
  )
} else if (capabilities("cairo")) {
  cairo_pdf(
    filename = output_file,
    width = plot_width_in,
    height = plot_height_in,
    family = "Arial"
  )
} else {
  pdf(
    file = output_file,
    width = plot_width_in,
    height = plot_height_in,
    family = "sans",
    useDingbats = FALSE
  )
}

print(p)
dev.off()
