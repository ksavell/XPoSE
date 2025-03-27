# F5
library(Seurat)
library(tidyverse)
library(venn)


source("Scripts/Functions/single_factor_DESeq.R")
source("Scripts/Functions/de_and_summary.R")
load("~/Projects/XPoSE/all_10312024.RData")

clusters <- unique(all$cluster_name)
all$experience <- ifelse(all$group == "Homecage", "HC", "NC")

pairs_list <- list(
  c("experience","NC","HC"),
  c("group", "Non-active", "Homecage"),
  c("group", "Active", "Homecage"),
  c("group", "Active", "Non-active")
)

de_and_summary(all,
               pair = pairs_list[[1]]) # this goes in Figure S4

de_and_summary(all,
               pair = pairs_list[[2]])

de_and_summary(all,
               pair = pairs_list[[3]])

de_and_summary(all,
               pair = pairs_list[[4]])


# compare AvHC and AvNA ---------------------------------------------------

source("Scripts/Functions/compare_deg_directories.R")

results <- compare_deg_directories(
  dir1 = "~/Projects/XPoSE/group_Active_Homecage",
  dir2 = "~/Projects/XPoSE/group_Active_Non-active",
  dir1_name = "ActivevHomecage",
  dir2_name = "ActivevNonactive",
  prefixes = clusters,  # Manually specify prefixes
  output_dir = "~/Projects/XPoSE/correspondence_AvHC_AvNA"
)

# deg overlap and venn ----------------------------------------------------

source("Scripts/Functions/deg_overlap_and_venn.R")

# Active v homecage
deg_overlap_and_venn(base_dir = "group_Active_Homecage", pattern = "Active_Homecage")

# Active v Nonactive
deg_overlap_and_venn(base_dir = "group_Active_Non-active", pattern = "Active_Non-active")


# DEG examples ------------------------------------------------------------

source("Scripts/Functions/summarize_de_counts.R")
source("Scripts/Functions/find_DEcounts.R")

gene_list <- c("Vgf","Ptprn", "Scg2", # common across most
               "Egr3","Fosb","Homer1","Nptx2", # all glut
               "Egr2", "Penk", "Mapk4", # all IT gluts
               "Arc","Bdnf", # notable arg's
               "Crh", #Pvalb specific
               "Cpne7","Ecel1" #Sst specific
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

