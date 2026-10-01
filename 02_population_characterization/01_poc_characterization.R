# Characterization for POC dataset

# Loading -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

source('02_population_characterization/functions/calc_prop.R')
source('02_population_characterization/functions/make_stdf.R')

# Paths -------------------------------------------------------------------------------
input_file <- 'output/01_metadata_clustering_qc/poc_hc_annotated.rds'
poc_hc <- readRDS(input_file)

output_dir <- 'output/02_population_characterization/poc'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Settings ----------------------------------------------------------------------------
cluster_colors <- c(
  'ITL23'      = '#2EBF5E',
  'ITL5'       = '#50B2AD',
  'ITL6'       = '#58D2CF',
  'ITvm'       = '#B1DE7D',
  'CTL6'       = '#2D8CB8',
  'CTL6b'      = '#7044AA',
  'ETL5'       = '#0D5A8B',
  'NPL5'       = '#3E9E64',
  'Pvalb'      = '#B9342C',
  'Sst'        = '#FF9900',
  'PvalbChand' = '#FF2D4E',
  'SstChodl'   = '#B1B10C',
  'Vip'        = '#B864CC',
  'Lamp5'      = '#DA808C',
  'Sncg'       = '#D3408D'
)

cluster_labels <- c(
  'ITL23'      = 'IT L2/3',
  'ITL5'       = 'IT L5',
  'ITL6'       = 'IT L6',
  'ITvm'       = 'IT vm',
  'CTL6'       = 'CT L6',
  'CTL6b'      = 'CT L6b',
  'ETL5'       = 'ET L5',
  'NPL5'       = 'NP L5',
  'Pvalb'      = 'Pvalb',
  'Sst'        = 'Sst',
  'PvalbChand' = 'Pvalb Chand',
  'SstChodl'   = 'Sst Chodl',
  'Vip'        = 'Vip',
  'Lamp5'      = 'Lamp5',
  'Sncg'       = 'Sncg'
)

cluster_order <- c(
  'ITL23', 
  'ITL5', 
  'ITL6', 
  'ITvm', 
  'CTL6', 
  'CTL6b', 
  'ETL5', 
  'NPL5', 
  'Pvalb', 
  'Sst', 
  'PvalbChand', 
  'SstChodl', 
  'Vip', 
  'Lamp5', 
  'Sncg'
)

# Cell type verification --------------------------------------------------------------
poc_hc_celltype <- calc_prop(seur_obj = poc_hc, 
                          fact1 = 'ratID',
                          fact2 = 'celltype',
                          fact3 = 'capture')

write.csv(poc_hc_celltype, file.path(output_dir, 'F1H_celltype_by_rat_capture.csv'))

# Mean reads by XPoSE-tag -------------------------------------------------------------
df <- make_stdf(poc_hc)

# Splitting the dataframe by both 'st' and 'cart'
IDslist <- split(df, list(df$st, df$cart))

# Calculating the mean for each split population
mean <- sapply(IDslist, function(x) {
  numeric_cols <- sapply(x, is.numeric)
  numeric_cols <- numeric_cols & !sapply(x, is.character)
  means <- apply(x[, numeric_cols], 2, mean, na.rm = TRUE)
  names(means) <- names(x)[numeric_cols]
  return(means)
})

write.csv(mean, file.path(output_dir, 'F1I_stReadsMean_bycart.csv'))

# Stats 
wilcox_results <- list()

# Explicit mapping between metadata tag labels and read-count columns
tag_map <- list(
  SampleTag08_mm = list(
    correct = "xpose_tag_08_reads",
    incorrect = c(
      "xpose_tag_02_reads",
      "xpose_tag_04_reads",
      "xpose_tag_06_reads"
    )
  ),
  SampleTag04_mm = list(
    correct = "xpose_tag_04_reads",
    incorrect = c(
      "xpose_tag_02_reads",
      "xpose_tag_06_reads",
      "xpose_tag_08_reads"
    )
  ),
  SampleTag02_mm = list(
    correct = "xpose_tag_02_reads",
    incorrect = c(
      "xpose_tag_04_reads",
      "xpose_tag_06_reads",
      "xpose_tag_08_reads"
    )
  ),
  SampleTag06_mm = list(
    correct = "xpose_tag_06_reads",
    incorrect = c(
      "xpose_tag_02_reads",
      "xpose_tag_04_reads",
      "xpose_tag_08_reads"
    )
  )
)

for (assigned_tag in names(tag_map)) {
  correct_reads_var <- tag_map[[assigned_tag]]$correct
  incorrect_reads_vars <- tag_map[[assigned_tag]]$incorrect

  # Cells assigned to this raw SampleTag label
  cells <- WhichCells(
    poc_hc,
    expression = xpose_tag == assigned_tag
  )

  # Correct-tag reads
  correct <- FetchData(
    poc_hc,
    vars = correct_reads_var
  )[cells, 1]

  # Mean reads across the three incorrect tags
  incorrect <- rowMeans(
    FetchData(
      poc_hc,
      vars = incorrect_reads_vars
    )[cells, , drop = FALSE]
  )

  # Paired Wilcoxon test
  test <- wilcox.test(
    correct,
    incorrect,
    paired = TRUE,
    alternative = "greater"
  )

  wilcox_results[[assigned_tag]] <- data.frame(
    sample = assigned_tag,
    correct_read_column = correct_reads_var,
    n_cells = length(cells),
    p_value = test$p.value,
    statistic = unname(test$statistic),
    method = test$method,
    alternative = test$alternative
  )
}

wilcox_summary <- do.call(rbind, wilcox_results)
write.csv(wilcox_summary, file.path(output_dir, "F1I_stats.csv"), row.names = FALSE)

# xpose_tag contribution / bias score per cluster ------------------------------------
# Pull metadata 
md <- poc_hc@meta.data %>%
  dplyr::select(cluster_name, xpose_tag, capture) %>%
  dplyr::mutate(
    cluster_name = as.character(cluster_name),
    xpose_tag   = as.character(xpose_tag),
    capture   = as.character(capture)
  )
tags     <- sort(unique(md$xpose_tag))
n_tags   <- length(tags)
expected <- 1 / n_tags          # equal-contribution null: 1/n_tags

# POOLED: Calculate counts + proportions per cluster 
counts <- md %>%
  dplyr::count(cluster_name, xpose_tag, name = 'n') %>%
  tidyr::complete(cluster_name, xpose_tag, fill = list(n = 0)) %>%
  dplyr::group_by(cluster_name) %>%
  dplyr::mutate(cluster_total = sum(n), prop = n / cluster_total) %>%
  dplyr::ungroup()

# Plots
counts$cluster_name      <- factor(counts$cluster_name,      levels = cluster_order)

# Formatted colored labels
axis_labels <- setNames(paste0("<span style='color:", cluster_colors[cluster_order], "; '>", 
                               cluster_labels[cluster_order], "</span>"), cluster_order)

# Pooled stacked bar
p_stack <- ggplot(counts, aes(x = prop, y = cluster_name, fill = xpose_tag)) +
  geom_col(width = 0.72, color = 'white', linewidth = 0.75) +
  scale_fill_grey(start = 0, end = 0.7) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_y_discrete(limits = rev(cluster_order), labels = axis_labels) +
  labs(x = 'Sample composition', y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = ggtext::element_markdown(size = 20, margin = margin(r = 10)),
    axis.text.x = element_text(size = 20, color = 'black'),
    axis.title.x = element_text(size = 24, margin = margin(t = 14)),
    axis.line = element_line(color = 'black', linewidth = 1),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(0.3, 'cm'),
    legend.position = 'none',
    plot.margin = margin(t = 15, r = 50, b = 15, l = 50))
ggsave(file.path(output_dir, 'stacked_bar_sampletag.pdf'),
       p_stack, width = 7, height = 5)
