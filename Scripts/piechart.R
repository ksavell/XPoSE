library(ggplot2)
library(dplyr)

hex_list <- c(
  'ITL23' = '#2EBF5E',
  'ITL5' = '#50B2AD',  
  'ITL6' = '#58D2CF',
  'ITvm' = '#B1DE7D',
  'CTL6' = '#2D8CB8',
  'CTL6b' = '#7044AA',
  'ETL5' = '#0D5A8B',
  'NPL5' = '#3E9E64',
  'Pvalb' = '#B9342C',
  'Sst' = '#FF9900',
  'PvalbChand' = '#FF2D4E',
  'SstChodl' = '#B1B10C',
  'Vip' = '#B864CC',
  'Lamp5' = '#DA808C',
  'Sncg' = '#D3408D'
)

all <- read.csv("all_cluster_proportions_long_07202026.csv")

all <- all %>%
  filter(region_analysis == "vmPFC")

axis_order <- c("ITL23", "ITL5", "ITL6", "ITvm", "CTL6", "CTL6b", "ETL5", "NPL5", "Pvalb", "Sst", "PvalbChand", "SstChodl", "Vip", "Lamp5", "Sncg")

experience_levels <- c(
  "N", "NT", "NC - Active", "NC - Non-active",
  "RT - Active", "RT - Non-active"
)

cluster_avg <- all %>%
  group_by(cluster_name, experience) %>%
  summarise(
    mean_percent = mean(percent * 100, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cluster_name = factor(cluster_name, levels = axis_order)
  )

plot_experiences <- intersect(
  experience_levels,
  unique(as.character(cluster_avg$experience))
)

for (exp_i in plot_experiences) {
  
  plot_df <- cluster_avg %>%
    filter(experience == exp_i) %>%
    arrange(cluster_name)
  
  pie_plot <- ggplot(
    plot_df,
    aes(x = "", y = mean_percent, fill = cluster_name)
  ) +
    geom_col(
      width = 1,
      color = "white",
      linewidth = 0.5,
      position = position_stack(reverse = FALSE)
    ) +
    coord_polar(
      theta = "y",
      direction = -1
    ) +
    scale_fill_manual(
      values = hex_list,
      breaks = axis_order,
      limits = axis_order,
      drop = FALSE
    ) +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(
    filename = paste0(
      gsub("[^A-Za-z0-9]+", "_", exp_i),
      "_vmPFC_pie_plot.pdf"
    ),
    plot = pie_plot,
    width = 6,
    height = 6,
    dpi = 300
  )
}