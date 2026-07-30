# Figure 5

# Loading -----------------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)
library()

source("Scripts/Functions/single_factor_DESeq.R")
source("Scripts/Functions/de_and_summary.R")
source("Scripts/Functions/compare_deg_directories.R")
source("Scripts/Functions/deg_overlap_and_venn.R")
source("Scripts/Functions/summarize_de_counts.R")
source("Scripts/Functions/find_DEcounts.R")

load("~/Projects/XPoSE/all_10312024.RData")


# F5A ---------------------------------------------------------------------

clusters <- unique(all$cluster_name)

pairs_list <- list(
  c("group", "Non-active", "Homecage"),
  c("group", "Active", "Homecage"),
  c("group", "Active", "Non-active")
)

de_and_summary(all,
               pair = pairs_list[[1]]) # NC:NA vs HC

de_and_summary(all,
               pair = pairs_list[[2]]) # NC:A vs HC

de_and_summary(all,
               pair = pairs_list[[3]]) # NC:A vs NC:NA


# F5B ---------------------------------------------------------------------

results <- compare_deg_directories(
  dir1 = "~/Projects/XPoSE/group_Active_Homecage",
  dir2 = "~/Projects/XPoSE/group_Active_Non-active",
  dir1_name = "ActivevHomecage",
  dir2_name = "ActivevNonactive",
  prefixes = clusters,  # Manually specify prefixes
  output_dir = "F5B_correspondence_AvHC_AvNA"
)


# F5C-D -------------------------------------------------------------------

deg_overlap_and_venn(base_dir = "group_Active_Homecage", pattern = "Active_Homecage")


# F5E ---------------------------------------------------------------------

clusters <- c("ITL23","ITL5","ITL6","CTL6","PTL5","Pvalb","Sst")

gene_list <- c("Vgf","Ptprn", "Scg2", # common across most
               "Egr3","Fosb","Homer1","Nptx2", # all glut
               "Penk", "Mapk4","Egr2", # all IT gluts
               "Arc","Bdnf", # notable arg's
               "Crh", #Pvalb specific
               "Ecel1" #Sst specific
)

for (cl in clusters) {
  find_DEcounts(directory = "group_Active_Homecage", 
                cluster = cl, 
                de_path = "_group_Active_Homecage", 
                feature_list = gene_list, 
                control_suffix = "HC")
}

#summarize for easy input into prism
summarize_DEcounts("group_Active_Homecage",
                   clusters = clusters,
                   gene_list = gene_list)


# DE updated --------------------------------------------------------------
summary_df <- as.data.frame(read_csv("HPC/output_local/output/downsample_vmPFC_07292026/de_summary_long_07292026.csv"))


cluster_cols = c('ITL23' = '#2EBF5E',
                 'ITL5' =	'#50B2AD',
                 'ITL6' =	'#58D2CF',
                 'CTL6' =	'#2D8CB8',
                 'ETL5' =	'#0D5A8B',
                 'Sst' =	'#FF9900'
                 )

cl_order <- c("ITL23", "ITL5", "ITL6", "CTL6", "ETL5", 
              "Sst") 
date_tag <- format(Sys.Date(), "%m%d%Y")

summary_df$cluster <- factor(summary_df$cluster, levels = cl_order)

p_plot <- ggplot(
  summary_df,
  aes(x = percentage, y = mean_DE, colour = cluster, group = cluster)
) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_DE - sd_DE), ymax = mean_DE + sd_DE),
    width = 3, linewidth = 0.25, alpha = 0.75
  ) +
  geom_line(linewidth = 0.25, alpha = 0.4) +
  geom_point(size = 0.25,
             shape = 16, alpha = 0.75) +
  scale_colour_manual(values = cluster_cols, breaks = cl_order, drop = FALSE) +
  scale_x_continuous(
    breaks = c(0, 1, 5, 25, 50, 75, 100),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", linewidth = 0.25) +
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
    legend.position = "bottom", # change to bottom etc to bring back
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
  file = paste0("de_by_active_fraction_", date_tag, ".pdf"),
  width = 2,
  height = 1.5,
  family = "Arial"
)
print(p_plot)
dev.off()

