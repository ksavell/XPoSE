# Cluster object

# Info --------------------------------------------------------------------

# This script creates
#       * glut subset object
#       * gaba subset object

# Loading -----------------------------------------------------------------
## Load packages -----------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)

## Data loading ------------------------------------------------------------
# load in initial 'combined' object that is output of create_object01.R

load("~/MakeObject/combined08032023.RData")

# First clustering ------------------------------------------------------

source("~/XPoSE/Scripts/Functions/cluster_first.R")

combined <- cluster_first(combined)

# Explore clusters --------------------------------------------------------

# First check QC metrics to see if any cluster is defined by low QC measures

source("~/XPoSE/Scripts/Functions/clstr_vln.R")

clstr_vln(combined, qc = T)

# Plot excitatory, inhibitory, and glia contamination markers
clstr_vln(combined, all = T)

# 26, 28, 30 are minor glia contamination
# 24 is low quality

# Subset glut -------------------------------------------------------------

#Specify which clusters to keep
keep_glut <- as.character(c(0:8,12:18,20,22,25))

source("~/XPoSE/Scripts/Functions/subset_reclust.R")

glut <- subset_reclust(combined, clust_tokeep = keep_glut, neigh_dim = 1:30, 
                       umap_dim = 1:30, res = 0.2)

# check glut clusters for known markers
clstr_vln(glut, glut = T)

# name the final clusters and save
glut_ids <- c("IT_L5/6", "IT_L2/3", "CT_L6", "PT_L5", "NP_L5/6", "CT_L6b")
names(glut_ids) <- levels(glut)
glut <- RenameIdents(glut, glut_ids)
glut$cluster_name <- paste(glut@active.ident)
save(glut, file = "glut.RData")


# Subset gaba -------------------------------------------------------------


keep_gaba <- as.character(c(9:11,19,21,23,27,32))

gaba <- subset_reclust(combined, clust_tokeep = keep_gaba, neigh_dim = 1:30, 
                       umap_dim = 1:30, res = 0.2)

# one more refinement

keep_gaba2 <- as.character(c(0:2,4:7))

gaba <- subset_reclust(gaba, clust_tokeep = keep_gaba2, neigh_dim = 1:25, 
                       umap_dim = 1:25, res = 0.2)

# check gaba clusters for known markers
clstr_vln(gaba, gaba = T)

# name final clusters and save
gaba_ids <- c("Pvalb", "Sst", "Ppp1r1b", "Vip", "Meis2", "Lamp5", "Sst Chodl")
names(gaba_ids) <- levels(gaba)
gaba <- RenameIdents(gaba, gaba_ids)
gaba$cluster_name <- paste(gaba@active.ident)
save(gaba, file = "gaba.RData")
