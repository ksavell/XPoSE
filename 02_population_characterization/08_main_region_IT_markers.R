# IT cell-type marker dot plot (vmPFC naive subset)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# Load in clustered Main object that is output of createobject_01.R
load('dmvmPFC_annotated_07162026.RData')

# User settings -----------------------------------------------------------------------

cluster_col <- 'cluster_name'
assay_use   <- 'RNA'
it_order <- c('ITvm', 'ITL6', 'ITL5', 'ITL23')

markers_by_cluster <- list(
  ITL23 = c('Enpp2', 'Glis3', 'Matn2', 'Cux2', 'Otof'),
  ITL5  = c('Rorb', 'Rmst', 'Zmat4', 'Scube1', 'Kcnab1'),
  ITL6  = c('Sema3e', 'Dcc', 'Galnt14', 'Sema3c', 'Col6a2'),
  ITvm  = c('Ndst4', 'Ndnf', 'Gipr', 'Trpc7', 'Pmp22', 'Tpbg')
)

expression_low    <- '#444444'
expression_mid    <- '#FFFFFF'
expression_high   <- '#000080'
expression_limits <- c(-2.5, 2.5)
expression_breaks <- c(-2.5, 0, 2.5)

output_root <- 'output'
save_png <- TRUE

pt <- function(x) x / 2.835

theme_pub <- theme_classic(base_size = 7, base_family = 'Arial') +
  theme(
    axis.line = element_line(linewidth = pt(0.5), colour = 'black'),
    axis.ticks = element_line(linewidth = pt(0.5), colour = 'black'),
    axis.text = element_text(size = 7, colour = 'black'),
    axis.title = element_text(size = 8, colour = 'black'),
    legend.position = 'right',
    legend.box = 'vertical',
    legend.title = element_text(size = 7, colour = 'black'),
    legend.text = element_text(size = 7, colour = 'black'),
    legend.key = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

# Plot data ---------------------------------------------------------------------------

# Subset naive vmPFC 
naive <- subset(all, subset = experience == 'N')
vm <- subset(naive, subset = region == 'vmPFC')

DefaultAssay(vm) <- assay_use
it <- subset(vm, subset = !!sym(cluster_col) %in% it_order)
it[[cluster_col]] <- factor(
  it[[cluster_col, drop = TRUE]],
  levels = it_order
)
Idents(it) <- it[[cluster_col, drop = TRUE]]

cat('cells per IT cluster (vmPFC naive):\n')
print(table(Idents(it)))

# Ordered gene list + presence check
genes <- unlist(markers_by_cluster, use.names = FALSE)
present <- genes %in% rownames(it[[assay_use]])

if (any(!present)) {
  cat(
    '\n** WARNING: genes not found and dropped:',
    paste(genes[!present], collapse = ', '),
    '**\n'
  )
}

genes_plot <- genes[present]
if (length(genes_plot) == 0) {
  stop('None of the requested marker genes were found.')
}

# Obtain Seurat DotPlot summaries, then rebuild manually 
dot_seed <- DotPlot(
  object = it,
  features = genes_plot,
  assay = assay_use,
  dot.scale = 5,
  col.min = expression_limits[1],
  col.max = expression_limits[2]
)

dot_data <- dot_seed$data %>%
  mutate(
    id = factor(as.character(id), levels = it_order),
    features.plot = factor(
      as.character(features.plot),
      levels = genes_plot
    ),
    avg.exp.scaled = pmax(
      pmin(avg.exp.scaled, expression_limits[2]),
      expression_limits[1]
    )
  )

# Final plot
p <- ggplot(
  dot_data,
  aes(x = features.plot, y = id)
) +
  geom_point(
    aes(size = pct.exp, fill = avg.exp.scaled),
    shape = 21,
    colour = 'grey35',
    stroke = pt(0.25)
  ) +
  scale_size_continuous(
    range = c(0, 5),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75),
    name = 'Percent expressed',
    guide = guide_legend(
      reverse = FALSE,
      order = 2,
      title.position = 'top'
    )
  ) +
  scale_fill_gradient2(
    low = expression_low,
    mid = expression_mid,
    high = expression_high,
    midpoint = 0,
    limits = expression_limits,
    breaks = expression_breaks,
    oob = scales::squish,
    na.value = 'grey90',
    name = 'Average\nexpression',
    guide = guide_colorbar(
      display = 'rectangles',
      reverse = FALSE,
      direction = 'vertical',
      order = 1,
      title.position = 'top',
      barheight = grid::unit(16, 'mm'),
      barwidth = grid::unit(2.5, 'mm')
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = 'italic'
    )
  )

# Save Quartz vector PDF + PNG --------------------------------------------------------

today <- format(Sys.Date(), '%m%d%Y')
out_dir <- file.path(output_root, paste0('IT_marker_dotplot_', today))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n <- length(genes_plot)
plot_width <- 4
plot_height <- 1.5

pdf_file <- file.path(
  out_dir,
  paste0('IT_marker_dotplot_', today, '.pdf')
)
png_file <- file.path(
  out_dir,
  paste0('IT_marker_dotplot_', today, '.png')
)

if (capabilities('aqua')) {
  quartz(
    type = 'pdf',
    file = pdf_file,
    width = plot_width,
    height = plot_height,
    family = 'Arial'
  )
} else {
  cairo_pdf(
    filename = pdf_file,
    width = plot_width,
    height = plot_height,
    family = 'Arial'
  )
}
print(p)
dev.off()

if (isTRUE(save_png)) {
  ggsave(
    filename = png_file,
    plot = p,
    width = plot_width,
    height = plot_height,
    dpi = 600,
    bg = 'white'
  )
}

write.csv(
  dot_data,
  file.path(
    out_dir,
    paste0('IT_marker_dotplot_plot_data_', today, '.csv')
  ),
  row.names = FALSE
)

