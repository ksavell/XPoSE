# Figure 5

# Loading -----------------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)

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
  output_dir = "~/Projects/XPoSE/correspondence_AvHC_AvNA"
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