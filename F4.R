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
