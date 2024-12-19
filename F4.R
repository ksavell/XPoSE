# Figure 4

# Info --------------------------------------------------------------------

# This script 


# Load libraries and data -------------------------------------------------

library(Seurat)
library(tidyverse)

load("all_10312024.RData")

source("Scripts/Functions/save_dimplot.R")

all$experience <- ifelse(all$group == "Homecage", "HC",
                          "NC")

# Custom groupby-hex combinations
hex_list <- list(
  'cluster_name' = c('CTL6' = '#2C8CB9',
                     'PTL5' = '#0A5B8C',
                     'ITL23' = '#41B75F',
                     'ITL5' = '#5DBFC1',  
                     'ITL6' = '#3A8F87',
                     'NPL56' = '#3C9E64',
                     'CTL6b' = '#6F499D',
                     'Pvalb' = '#E66027',
                     'Sst' = '#F8991D',
                     'Meis2' = '#C52126',
                     'Vip' = '#A669AB',
                     'Lamp5' = '#DB808C',
                     'SstChodl' = '#B0B235',
                     'PvalbChand' = '#AD6C49'),
  'orig.ident' = c('C1' = '#22677A', 
                'C2' = '#D4A841'),
  'sex' = c('male' = "#2C5F2D", 
            'female' = "#97BC62"),
  'experience' =  c('HC' = '#c0c0c0',
                    'NC' = '#ae1e5b'),
  'group' = c('Non-active' = '#e37a9e',
              'Active' = '#801743')
)

save_dimplot(all, 
             groupby = 'cluster_name',
             file_n = 'all',
             hex_list = hex_list)

save_dimplot(all, 
             groupby = 'experience',
             file_n = 'all',
             hex_list = hex_list)

# nc plots 

nc <- subset(all, subset = experience == 'NC')

save_dimplot(nc, 
             groupby = 'group',
             splitby = 'ratID',
             file_n = 'nc',
             hex_list = hex_list)

# Cluster prop by rat/capture ---------------------------------------------

source("Scripts/Functions/calc_prop.R")

obj_celltype <- calc_prop(seur_obj = all, 
                          fact1 = 'ratID',
                          fact2 = 'group')

write.csv(obj_celltype, file = "all_counts_by_rat_group.csv")

clust_prop_cart <- calc_prop(seur_obj = all, 
                             fact1 = 'ratID',
                             fact2 = 'cluster_name',
                             fact3 = 'group')

write.csv(clust_prop_cart, file = "all_clust_prop_group.csv")



# now plot the data -------------------------------------------------------

plot_data <- all_group_cluster %>%
  group_by(cluster_name, group) %>%
  summarize(mean_percent = mean(percent, na.rm = TRUE),
            sd_percent = sd(percent, na.rm = TRUE),
            mean_count = mean(count, na.rm = TRUE),
            sd_count = sd(count, na.rm = TRUE),
            .groups = "drop")

all_group_cluster$group <- factor(all_group_cluster$group, levels = c("Homecage", "Non-active", "Active"))
plot_data$group <- factor(plot_data$group, levels = c("Homecage", "Non-active", "Active"))

# Define custom colors for each group
group_colors <- c("Homecage" = "#C0C0C0", "Non-active" = "#e37a9e", "Active" = "#801743")
cluster_colors <- c('CTL6' = '#2C8CB9',
                    'PTL5' = '#0A5B8C',
                    'ITL23' = '#41B75F',
                    'ITL5' = '#5DBFC1',  
                    'ITL6' = '#3A8F87',
                    'NPL56' = '#3C9E64',
                    'CTL6b' = '#6F499D',
                    'Pvalb' = '#E66027',
                    'Sst' = '#F8991D',
                    'Meis2' = '#C52126',
                    'Vip' = '#A669AB',
                    'Lamp5' = '#DB808C',
                    'SstChodl' = '#B0B235',
                    'PvalbChand' = '#AD6C49')


# Plot with ggplot2
pdf("F4_prop_cluster_group.pdf",
    width = 6,
    height = 1.5)
ggplot(plot_data, aes(x = cluster_name, y = mean_percent, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_percent, ymax = mean_percent + sd_percent), 
                position = position_dodge(width = 0.5), width = 0, color = "black") +
  geom_jitter(data = all_group_cluster, aes(x = cluster_name, y = percent, color = group),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), size = 1, shape = 1, color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(y = "Proportion\n(% n[cluster] / n[group])",
       x = NULL)+
  theme_classic() +
  theme(
    axis.title.y = element_text(angle = 90, size = 8),
    axis.text.x = element_text(angle = 0, size = 8, color = cluster_colors),
    legend.position = "none"
  )
dev.off()
 

pdf("F4_count_cluster_group.pdf",
    width = 6,
    height = 1.5)
ggplot(plot_data, aes(x = cluster_name, y = mean_count, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_count, ymax = mean_count + sd_count), 
                position = position_dodge(width = 0.75), width = 0, color = "black") +
  geom_jitter(data = all_group_cluster, aes(x = cluster_name, y = count, color = group),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1, shape = 1, color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(y = "Proportion\n(% n[cluster] / n[group])",
       x = NULL)+
  theme_classic() +
  theme(
    axis.title.y = element_text(angle = 90,size = 8),
    axis.text.x = element_text(angle = 0, size = 8, color = cluster_colors),
    legend.position = "none"
  )
dev.off()



# calculate z score -------------------------------------------------------

library(dplyr)

calculate_z_scores_per_cluster <- function(data, numeric_column, group_column, celltype_column, control_group_name) {
  # Initialize an empty dataframe to store results
  results_df <- data.frame()
  
  # Get unique combinations of group and cell type
  clusters <- unique(data[c(group_column, celltype_column)])
  
  # Loop through each cluster
  for (i in seq(nrow(clusters))) {
    cluster_info <- clusters[i, ]
    cluster_data <- data[data[[group_column]] == cluster_info[[group_column]] & data[[celltype_column]] == cluster_info[[celltype_column]], ]
    
    if (cluster_info[[group_column]] == control_group_name) {
      # Calculate mean and SD for the control group within this cluster
      control_mean <- mean(cluster_data[[numeric_column]], na.rm = TRUE)
      control_sd <- sd(cluster_data[[numeric_column]], na.rm = TRUE)
    } else {
      # Assume control data is derived from the same cell type within the Homecage group
      control_data <- data[data[[group_column]] == control_group_name & data[[celltype_column]] == cluster_info[[celltype_column]], ]
      control_mean <- mean(control_data[[numeric_column]], na.rm = TRUE)
      control_sd <- sd(control_data[[numeric_column]], na.rm = TRUE)
    }
    
    # Calculate z-scores if valid control data exists
    if (!is.na(control_mean) && control_sd > 0) {
      cluster_data$z_score <- (cluster_data[[numeric_column]] - control_mean) / control_sd
      results_df <- rbind(results_df, cluster_data)
    }
  }
  
  return(results_df)
}


group_cluster_zscore <- calculate_z_scores_per_cluster(clust_prop_cart,
                                                       numeric_column = "count",
                                                       group_column = "group",
                                                       celltype_column = "cluster_name",
                                                       control_group_name = "Homecage")

group_means <- group_cluster_zscore %>%
  group_by(group, cluster_name) %>%
  summarize(mean_z_score = mean(z_score, na.rm = TRUE), .groups = 'drop')

# Create the plot
ggplot(group_means, aes(x = cluster_name, y = mean_z_score, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +  # Bar plot with dodging
  geom_point(data = group_cluster_zscore, aes(x = cluster_name, y = z_score, color = group),
             position = position_dodge(width = 0.8), size = 1.5, alpha = 0.6) +  # Add jittered points
  labs(title = "Average Z-Scores and Individual Data Points by Cluster and Group",
       x = "Cluster",
       y = "Z-Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better visibility
  scale_fill_brewer(palette = "Set1") +  # Color palette for bars
  scale_color_brewer(palette = "Set1")


# trying fold change angle now --------------------------------------------

library(dplyr)

calculate_fold_change_per_sample <- function(data, numeric_column, group_column, cluster_column, control_group_name) {
  # Check if columns exist in the dataframe
  if (!numeric_column %in% names(data) || !group_column %in% names(data) || !cluster_column %in% names(data)) {
    stop("One or more specified columns do not exist in the dataframe.")
  }
  
  # Initialize an empty dataframe to store results
  results_df <- data.frame()
  
  # Get unique clusters from the cluster column
  clusters <- unique(data[[cluster_column]])
  
  # Loop through each unique cluster
  for (cluster in clusters) {
    # Filter data for the current cluster
    cluster_data <- filter(data, !!sym(cluster_column) == cluster)
    
    # Calculate control group mean for the current cluster
    control_data <- filter(cluster_data, !!sym(group_column) == control_group_name)
    control_mean <- mean(control_data[[numeric_column]], na.rm = TRUE)
    
    # Check if control mean is valid
    if (is.na(control_mean) || control_mean == 0) {
      next  # Skip this cluster if control mean is NA or zero
    }
    
    # Calculate fold change for each sample in the cluster
    cluster_data$fold_change <- cluster_data[[numeric_column]] / control_mean
    
    # Append to results dataframe
    results_df <- rbind(results_df, cluster_data)
  }
  
  return(results_df)
}

fold_change_results <- calculate_fold_change_per_sample(clust_prop_cart, "percent", "group", "cluster_name", "Homecage")

write_csv(fold_change_results, "clust_prop_group_FC.csv")

# Assuming 'results_df' is your dataframe from the fold change calculation
# Calculate means for each group and cluster combination
group_means <- fold_change_results %>%
  group_by(group, cluster_name) %>%
  summarize(mean_fold_change = mean(fold_change, na.rm = TRUE), .groups = 'drop')

# Create the plot
ggplot(group_means, aes(x = cluster_name, y = mean_fold_change, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +  # Bar plot with dodging
  geom_point(data = fold_change_results, aes(x = cluster_name, y = fold_change, color = group),
             position = position_dodge(width = 0.8), size = 1.5, alpha = 0.6) +  # Add jittered points
  labs(title = "Average Fold Change and Individual Data Points by Cluster and Group",
       x = "Cluster",
       y = "Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better visibility
  scale_fill_brewer(palette = "Set1") +  # Color palette for bars
  scale_color_brewer(palette = "Set1")  # Color palette for points
