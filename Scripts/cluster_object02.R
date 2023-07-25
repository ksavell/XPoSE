# Cluster object

# By: Katherine Savell


# Info --------------------------------------------------------------------

# This script creates
#       * glut subset object
#       * gaba subset object


# Loading -----------------------------------------------------------------
## Load packages -----------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)
#library(here)

## Data loading ------------------------------------------------------------
# load in initial 'combined' object that is output of create_object01.R

load("~/MakeObject/combined07202023.RData")

# Remove pseudogenes ------------------------------------------------------

source("~/XPoSE/Scripts/Functions/remove_genes.R")

combined <- remove_genes(combined)

# First clustering ------------------------------------------------------

source("~/XPoSE/Scripts/Functions/cluster_first.R")

cluster_first(combined)

# Explore clusters --------------------------------------------------------

# First check QC metrics to see if any cluster is defined by low QC measures

source("~/XPoSE/Scripts/Functions/clstr_vln.R")

clstr_vln(combined, qc = T)

# Plot excitatory, inhibitory, and glia contamination markers
clstr_vln(combined, all = T)

# Subset glut and gaba ----------------------------------------------------

## In progress ##

#Specify which clusters to keep
keep_glut <- as.character(c(0:8,12:17,19:21,24,29))

source("~/XPoSE/Scripts/Functions/subset_reclust.R")

glut <- subset_reclust(combined, clust_tokeep = keep_glut, neigh_dim = 1:30, umap_dim = 1:30, res = 0.2)

#this is what this should be doing
glut <- subset(combined, idents = c(as.character(0:8), as.character(12:17), 
                                    as.character(19:21), "24", "29"))

glut <- glut %>%
        FindVariableFeatures(selection.method = "vst", 
                             nfeatures = 2000) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(ndims.print = 1:50,
               nfeatures.print = 50)

glut <- glut %>%
        FindNeighbors(dims = 1:30) %>%
        RunUMAP(reduction = "pca", 
                dims = 1:30)

glut <- FindClusters(glut, 
                         resolution = 0.2)

glut_ids <- c("IT_L5/6", "CT_L6", "IT_L2/3", "PT_L5", "NP_L5/6", "IT_L2/3", "CT_L6b")
names(glut_ids) <- levels(glut)
glut <- RenameIdents(glut, glut_ids)

save(glut, file = "glut07202023.RData")

## gaba time

gaba <- subset(combined, idents = c(as.character(9:11), "18", as.character(22:23), 
                                    "26", as.character(31:32)))

gaba <- gaba %>%
        FindVariableFeatures(selection.method = "vst", 
                             nfeatures = 2000) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(ndims.print = 1:50, 
               nfeatures.print = 50)

gaba <- gaba %>%
        FindNeighbors(dims = 1:50) %>%
        RunUMAP(reduction = "pca", 
                dims = 1:50)

gaba <- FindClusters(gaba, 
                     resolution = 2)

# doublet clusters 6 and 16,  exclude then recluster
gaba <- subset(gaba, idents = c(as.character(0:5),
                                as.character(7:15),
                                as.character(17:18)))
# recluster 0.2 resolution

gaba <- gaba %>%
        FindVariableFeatures(selection.method = "vst", 
                             nfeatures = 2000) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(ndims.print = 1:50, 
               nfeatures.print = 50)

gaba <- gaba %>%
        FindNeighbors(dims = 1:50) %>%
        RunUMAP(reduction = "pca", 
                dims = 1:50)

gaba <- FindClusters(gaba, 
                     resolution = 0.2)


#VlnPlot(gaba, features = c("Gad1","Kcnc2","Sst","Chodl","Ppp1r1b", "Meis2","Vip","Lamp5", "Ndnf","Cck"))

gaba_ids <- c("Pvalb", "Sst", "Ppp1r1b", "Vip", "Meis2", "Lamp5", "Sst Chodl")
names(gaba_ids) <- levels(gaba)
gaba <- RenameIdents(gaba, gaba_ids)
save(gaba, file = "gaba07202023.RData")
