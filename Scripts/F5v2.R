# f5, group/rat downsample

source("Scripts/Functions/single_factor_DESeq.R")
load("~/Projects/XPoSE/all_10312024.RData")


clusters <- unique(all$cluster_name)
all$experience <- ifelse(all$group == "Homecage", "HC", "NC")

pairs_list <- list(
  c("experience","NC","HC"),
  c("group", "Non-active", "Homecage"),
  c("group", "Active", "Homecage"),
  c("group", "Active", "Non-active")
)

# Set seed
set.seed(22)

# Get metadata
meta <- all@meta.data

# Set the grouping factors
grouping <- meta %>%
  group_by(group, ratID) %>%
  summarise(n_cells = n(), .groups = "drop")

# Determine minimum cells across group-rat combos (372 in your case)
target_size <- min(grouping$n_cells)  # Or set manually: target_size <- 372

cat("Downsampling each group-ratID pair to", target_size, "cells.\n")

# Sample cells from each group-ratID combo
selected_cells <- grouping %>%
  rowwise() %>%
  mutate(cells = list({
    cells_in_group <- rownames(meta)[meta$group == group & meta$ratID == ratID]
    if (length(cells_in_group) >= target_size) {
      sample(cells_in_group, target_size)
    } else {
      # Optionally skip or include all if under target
      character(0)
    }
  })) %>%
  pull(cells) %>%
  unlist()

# Subset the Seurat object
all <- subset(all, cells = selected_cells)

# Check the final counts
table(all$group, all$ratID)


# overlap between AvHC and AvNA -------------------------------------------

AHC_dir <- ()

ANA_dir <- ()

