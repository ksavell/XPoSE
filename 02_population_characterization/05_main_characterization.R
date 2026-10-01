# Characterization for Main dataset

# Loading -----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)
library(scales)

source('02_population_characterization/functions/calc_prop.R')

# Paths -------------------------------------------------------------------------------
input_file <- 'output/02_population_characterization/main_annotated.rds'
main <- readRDS(input_file)

output_dir <- 'output/02_population_characterization/main'
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

# Counts by experience ----------------------------------------------------------------
obj_counts <- calc_prop(seur_obj = main, 
                          fact1 = 'ratID',
                          fact2 = 'experience')

obj_counts <- obj_counts[obj_counts$count != 0, ]

write.csv(obj_counts, file.path(output_dir, 'F3E_counts_by_experience.csv'))

# Activity defined captures -----------------------------------------------------------
experience_set <- subset(main, subset = experience == 'NC' | experience == 'RT')

activity_percent <- calc_prop(seur_obj = experience_set, 
                        fact1 = 'ratID',
                        fact2 = 'experience',
                        fact3 = 'population')

activity_percent <- activity_percent %>%
  filter(count != 0) %>%
  group_by(ratID, experience) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()

write.csv(activity_percent, file.path(output_dir, 'F3F_activity_defined_captures.csv'))

# XPoSE-tag composition across capture ------------------------------------------------
tag_composition <- calc_prop(seur_obj = main, 
                        fact1 = 'ratID',
                        fact2 = 'capture')

tag_composition <- tag_composition %>%
  filter(count != 0) %>%
  group_by(capture) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()

write.csv(tag_composition, file.path(output_dir, 'F3G_tag_composition_capture.csv'))

# Individual contribution / bias score per cluster ------------------------------------
# Subset by experience
exp_subset <- subset(x = main, subset = experience == 'NT')

# Pull metadata 
md <- exp_subset@meta.data %>%
  dplyr::select(cluster_name, ratID, capture) %>%
  dplyr::mutate(
    cluster_name = as.character(cluster_name),
    ratID   = as.character(ratID),
    capture   = as.character(capture)
  )
tags     <- sort(unique(md$ratID))
n_tags   <- length(tags)
expected <- 1 / n_tags          # equal-contribution null: 1/n_tags

# Calculate counts + proportions per cluster
counts <- md %>%
  dplyr::count(cluster_name, ratID, name = 'n') %>%
  tidyr::complete(cluster_name, ratID, fill = list(n = 0)) %>%
  dplyr::group_by(cluster_name) %>%
  dplyr::mutate(cluster_total = sum(n), prop = n / cluster_total) %>%
  dplyr::ungroup()

counts$cluster_name      <- factor(counts$cluster_name,      levels = cluster_order)

# Formatted colored labels
axis_labels <- setNames(paste0("<span style='color:", cluster_colors[cluster_order], "; '>", 
                               cluster_labels[cluster_order], "</span>"), cluster_order)

# Pooled stacked bar
p_stack <- ggplot(counts, aes(x = prop, y = cluster_name, fill = ratID)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.75) +
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
ggsave(file.path(output_dir, 'stacked_bar_rat.pdf'), p_stack, 
       width = 7, height = 5)
