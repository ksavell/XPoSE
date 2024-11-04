# Cluster object

# Info --------------------------------------------------------------------

# This script creates
#       * hc only object
#       * hc/nc combined object

# Loading -----------------------------------------------------------------
## Load packages -----------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)

## Data loading ------------------------------------------------------------
# load in initial 'combined' object that is output of create_object01.R

load("combined08032023.RData")

# Assign MapMyCells metadata ----------------------------------------------

mapping <- read.csv("combinedcounts_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1726072259014.csv",
                    comment.char = "#")
head(data.frame(mapping))

combined$class <- mapping$class_name
combined$subclass <- mapping$subclass_name

# ID neurons vs non-neurons

# Extract the class column from the metadata
metadata <- combined@meta.data

# Extract the numeric part of the 'class' column
metadata$numeric_part <- as.numeric(sub(" .*", "", metadata$class))

# Assign 'other' if the numeric part is 30 or greater, else 'neuron'
metadata$celltype <- ifelse(metadata$numeric_part >= 30, 'other', 'neuron')

# Update the Seurat object's metadata
combined@meta.data <- metadata

# Optionally, drop the helper column 'numeric_part'
combined@meta.data$numeric_part <- NULL

save(combined, file = "combined_withmmmannotation_09112024.RData")

# Subset out other cell types  ------------------------------

combined <- subset(combined, subset = celltype == 'neuron')

# Explore homecage clusters ----------------------------------------------

combined <- subset(combined, subset = group == 'Homecage')

source("Scripts/Functions/cluster_first.R")

combined <- cluster_first(combined)

# First check QC metrics to see if any cluster is defined by low QC measures

source("Scripts/Functions/clstr_vln.R")

clstr_vln(combined, qc = T)

# Plot excitatory, inhibitory, and glia contamination markers
clstr_vln(combined, all = T)

# 26 is minor glia contamination
# 21 is clustering by lower QC measures

keep <- as.character(c(0:20,22:25,27))

source("Scripts/Functions/subset_reclust.R")

combined_f <- subset_reclust(combined, clust_tokeep = keep, neigh_dim = 1:30, 
                       umap_dim = 1:30, res = 0.5)

keep2 <- as.character(c(0:12)) # got rid of Gja1 cluster

hc <- subset_reclust(combined_f, clust_tokeep = keep2, neigh_dim = 1:30, 
                             umap_dim = 1:30, res = 0.5)

hc_ids <- c("CTL6","ITL23", "ITL5","PTL5", "ITL6","Sst","Meis2","Pvalb",
            "NPL56","ITL23","CTL6b","Vip","Lamp5","SstChodl","PvalbChand")

names(hc_ids) <- levels(hc)
hc <- RenameIdents(hc, hc_ids)
hc$cluster_name <- paste(hc@active.ident)
save(hc, file = "hc_09112024.RData")


# Combined clustering -----------------------------------------------------
load("~/rstudio_projects/XPoSE/combined_withmmmannotation_09112024.RData")

combined <- subset(combined, subset = celltype == 'neuron')


source("Scripts/Functions/cluster_first.R")

combined <- cluster_first(combined)

# First check QC metrics to see if any cluster is defined by low QC measures

source("Scripts/Functions/clstr_vln.R")

clstr_vln(combined, qc = T)

# Plot excitatory, inhibitory, and glia contamination markers
clstr_vln(combined, all = T)

# 27 is minor glia contamination
# 28 is clustering by lower QC measures

keep <- as.character(c(0:26,29:31))

source("Scripts/Functions/subset_reclust.R")

combined_f <- subset_reclust(combined, clust_tokeep = keep, neigh_dim = 1:30, 
                             umap_dim = 1:30, res = 0.5)

clstr_vln(combined_f, qc = T)
clstr_vln(combined_f, all = T)

keep2 <- as.character(c(0:12,14:15)) # got rid of gaba mural cluster 13

all_f <- subset_reclust(combined_f, clust_tokeep = keep2, neigh_dim = 1:30, 
                     umap_dim = 1:25, res = 1) 

# get rid of one cluster of suspected duplicates
all <- subset_reclust(all_f, clust_tokeep = as.character(c(0:17,19)), neigh_dim = 1:40, 
                      umap_dim = 1:40, res = 1)

clstr_vln(all, qc = T)
clstr_vln(all, all = T)

all_ids <- rep(NA_character_, 19) # add 1 doesnt start at 0
all_ids[c(2,5,13)] <- "ITL23"
all_ids[c(1,15,17)] <- "ITL5"
all_ids[c(7)] <- "ITL6"
all_ids[c(6,11)] <- "PTL5"
all_ids[12] <- "NPL56"
all_ids[c(3,4)] <- "CTL6"
all_ids[14] <- "CTL6b"
all_ids[10] <- "Meis2"
all_ids[9] <- "Sst"
all_ids[8] <- "Pvalb"
all_ids[16] <- "Vip"
all_ids[18] <- "Lamp5"
all_ids[19] <- "SstChodl"
all_ids[20] <- "PvalbChand"

names(all_ids) <- levels(all)
all <- RenameIdents(all, all_ids)
all$cluster_name <- paste(all@active.ident)
save(all, file = "all_10312024.RData")

