library(tidyverse)
library(dplyr)

all <- read.csv("all_cluster_proportions_long_07202026.csv")

all <- all %>%
  filter(region_analysis == "dmPFC")

axis_order <- c("ITL23", "ITL5", "ITL6", "ITvm", "CTL6", "CTL6b", "ETL5", "NPL5", "Pvalb", "PvalbChand", "Sst", "SstChodl", "Sncg", "Vip", "Lamp5")

experience_levels <- c(
  "N", "NT", "NC - Active", "NC - Non-active",
  "RT - Active", "RT - Non-active"
)

# experience_levels <- c(
#   "HC", "NC - Active", "NC - Non-active"
# )

all <- all %>%
  mutate(
    cell_type = if_else(
      cluster_name %in% c(
        "ITL23", "ITL5", "ITL6", "ITvm",
        "CTL6", "CTL6b", "ETL5", "NPL5"
      ),
      "excitatory",
      "inhibitory"
    )
  )

cluster_colors <- c(
  'CTL6' = '#2D8CB8',
  'CTL6b' = '#7044AA',
  'ETL5' = '#0D5A8B',
  'ITL23' = '#2EBF5E',
  'ITL5' = '#50B2AD',  
  'ITL6' = '#58D2CF',
  'ITvm' = '#B1DE7D',
  'NPL5' = '#3E9E64',
  'Pvalb' = '#B9342C',
  'PvalbChand' = '#FF2D4E',
  'Sst' = '#FF9900',
  'SstChodl' = '#B1B10C',
  'Sncg' = '#D3408D',
  'Vip' = '#B864CC',
  'Lamp5' = '#DA808C'
)

tier_colors <- c(
  "10%" = "grey70",
  "20%" = "grey50",
  "30%" = "grey30"
    
)

experience_colors <- c(
  "N" = "black",
  "NT" = "grey70",
  "NC - Active" = "#00675F",
  "NC - Non-active" = "#75C3BC",
  "RT - Active" = "#801743",
  "RT - Non-active" = "#e37a9e"
)

# experience_colors <- c(
#   "HC" = "grey70",
#   "NC - Active" = "#00675F",
#   "NC - Non-active" = "#75C3BC"
# )

# Summarize
summary_df <- all %>%
  mutate(
    cluster = factor(cluster_name, levels = axis_order),
    experience = factor(experience, levels = experience_levels)
  ) %>%
  group_by(cell_type, experience, cluster) %>%
  summarise(
    mean = mean(proportion, na.rm = TRUE),
    sd = sd(proportion, na.rm = TRUE),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = pmax(mean - sem, 0),
    ymax = pmin(mean + sem, 1)
  )

cap_width <- 0.010
tier_levels <- c(0.10, 0.20, 0.30)
label_offset <- 0.01
plot_scale <- 1.5

plot_experiences <- intersect(
  experience_levels,
  unique(as.character(summary_df$experience))
)

plot_celltypes <- unique(summary_df$cell_type)

for (exp_i in plot_experiences) {
  for (cell_i in plot_celltypes) {
    
    exp_color <- experience_colors[[exp_i]]
  
    cell_axis_order <- axis_order[
      axis_order %in%
        as.character(
          summary_df$cluster[summary_df$cell_type == cell_i]
        )
    ]
    
    n_axes <- length(cell_axis_order)
    
    # Tier outlines for the current cell type
    guide_df <- map_dfr(tier_levels, function(r) {
      
      tier_name <- paste0(r * 100, "%")
      
      tibble(
        tier = tier_name,
        angle = 2 * pi *
          (seq_along(cell_axis_order) - 1) / n_axes,
        x = plot_scale * r * sin(angle),
        y = plot_scale * r * cos(angle),
        plot_color = unname(tier_colors[tier_name])
      ) %>%
        slice(c(seq_len(n()), 1))
    })
    
    # Spokes for the current cell type
    spoke_df <- tibble(
      cluster = factor(
        cell_axis_order,
        levels = cell_axis_order
      ),
      angle = 2 * pi *
        (seq_along(cell_axis_order) - 1) / n_axes,
      xend = plot_scale * 0.35 * sin(angle),
      yend = plot_scale * 0.35 * cos(angle),
      plot_color = unname(cluster_colors[cell_axis_order])
    )
    
    # Labels for the current cell type
    lab_df <- tibble(
      cluster = factor(cell_axis_order, levels = cell_axis_order),
      angle = 2 * pi * (seq_along(cell_axis_order) - 1) / n_axes,
      
      x = (0.523 + label_offset) * sin(angle),
      y = (0.523 + label_offset) * cos(angle),
      
      label = cell_axis_order,
      plot_color = unname(cluster_colors[cell_axis_order]),
      
       hjust = case_when(
    sin(angle) >  0.15 ~ 0,  
    sin(angle) < -0.15 ~ 1,  
    TRUE ~ 0.5                
  ),

  vjust = case_when(
    cos(angle) >  0.15 ~ 0, 
    cos(angle) < -0.15 ~ 1, 
    TRUE ~ 0.5                
  )
)
    
    # Plot data using the current cell-type axis order
    plot_df <- summary_df %>%
      filter(
        experience == exp_i,
        cell_type == cell_i,
        cluster %in% cell_axis_order
      ) %>%
      mutate(
        cluster = factor(
          as.character(cluster),
          levels = cell_axis_order
        )
      ) %>%
      arrange(cluster) %>%
      mutate(
        plot_color = unname(
          cluster_colors[as.character(cluster)]
        ),
        
        angle = 2 * pi *
          (as.numeric(cluster) - 1) / n_axes,
        
        x = plot_scale * mean * sin(angle),
        y = plot_scale * mean * cos(angle),
        
        x_min = plot_scale * ymin * sin(angle),
        y_min = plot_scale *  ymin * cos(angle),
        
        x_max = plot_scale * ymax * sin(angle),
        y_max = plot_scale * ymax * cos(angle),
        
        dx = plot_scale * cap_width * cos(angle),
        dy = plot_scale * -cap_width * sin(angle)
      )
    
    p <- ggplot() +
      
      geom_polygon(
        data = guide_df,
        aes(
          x = x,
          y = y,
          group = tier,
          color = plot_color,
          linewidth = tier
        ),
        fill = NA
      ) +
      
      geom_segment(
        data = spoke_df,
        aes(
          x = 0,
          y = 0,
          xend = xend,
          yend = yend,
          color = plot_color
        ),
        linetype = "dotted",
        linewidth = 0.6
      ) +
      
      geom_text(
        data = lab_df,
        aes(
          x = x,
          y = y,
          label = label,
          hjust = hjust,
          vjust = vjust,
          color = plot_color
        ),
        size = 10
      ) +
      
      geom_polygon(
        data = plot_df,
        aes(
          x = x,
          y = y,
          group = 1
        ),
        fill = adjustcolor(
          exp_color,
          alpha.f = 0.50
        ),
        color = exp_color,
        linewidth = 0
      ) +
      
      geom_segment(
        data = plot_df,
        aes(
          x = x_min,
          y = y_min,
          xend = x_max,
          yend = y_max
        ),
        color = "black",
        linewidth = 0.6
      ) +
      
      geom_segment(
        data = plot_df,
        aes(
          x = x_max - dx,
          y = y_max - dy,
          xend = x_max + dx,
          yend = y_max + dy
        ),
        color = "black",
        linewidth = 0.6
      ) +
      
      geom_segment(
        data = plot_df,
        aes(
          x = x_min - dx,
          y = y_min - dy,
          xend = x_min + dx,
          yend = y_min + dy
        ),
        color = "black",
        linewidth = 0.6
      ) +
      
      geom_point(
        data = plot_df,
        aes(
          x = x,
          y = y,
          color = plot_color
        ),
        shape = 21,
        fill = "white",
        stroke = 0.8,
        size = 5
      ) +
      
      scale_color_identity() +
      
      scale_linewidth_manual(
        values = c(
          "10%" = 0.45,
          "20%" = 0.60,
          "30%" = 0.75
        ),
        guide = "none"
      ) +
      
      coord_equal(
        xlim = c(-0.75, 0.75),
        ylim = c(-0.75, 0.75),
        clip = "off"
      ) +
      
      labs(
        title = paste(exp_i, cell_i, sep = " – ")
      ) +
      
      theme_void(base_size = 14) +
      
      theme(
        legend.position = "none",
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        ),
        plot.margin = margin(20, 100, 20, 40)
      )
    
    ggsave(
      filename = paste0(
        gsub("[^A-Za-z0-9]+", "_", exp_i),
        "_",
        cell_i,
        "_radial_plot.pdf"
      ),
      plot = p,
      width = 6,
      height = 6,
      dpi = 300
    )
  }
}