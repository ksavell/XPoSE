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
library(here)


## Data loading ------------------------------------------------------------
# load in the initial 'combined' object


# Functions ---------------------------------------------------------------
# Who knows what will end up here!
# Perhaps plot formatting? 
# remove pseudogenes?
# 

# Remove pseudogenes ------------------------------------------------------
# Get gene lists
genes <- combined@assays$RNA@counts@Dimnames[[1]]

# search known pseudogene, microRNA, model gene, or other non-coding patterns
indices_remove <- c(grep(pattern = "BC0", x = genes),
                    grep(pattern = ".ps", x = genes),
                    grep(pattern = "Rik", x = genes),
                    grep(pattern = "rik", x = genes),
                    grep(pattern = "Gm[0-9]", x = genes),
                    grep(pattern = "Mir[0-9]", x = genes),
                    grep(pattern = "LOC[0-9]", x = genes))

# remove genes
genes[-indices_remove] -> genes

#subset here again to save the day :') 
combined <- subset(combined, features = genes)

# should we re-calculate the nCount and nFeatures? 

# First clustering ------------------------------------------------------

combined <- combined %>%
        NormalizeData() %>%
        FindVariableFeatures(selection.method = "vst", 
                             nfeatures = 2000) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(ndims.print = 1:50, 
               nfeatures.print = 50) # printed features that make up PC, 
                                     # can take this out

# 50 PCAs actually seems appropriate, 
# I noticed gaba markers subtleties even in the 40s
combined <- combined %>%
        FindNeighbors(dims = 1:50) %>%
        RunUMAP(reduction = "pca", 
                dims = 1:50)

#starting with a high resolution
combined <- FindClusters(combined, 
                         resolution = 2)

# Explore clusters --------------------------------------------------------
# First check QC metrics to see if any cluster is defined by low QC measures
VlnPlot(combined, features = c("nFeature_RNA","nCount_RNA"))

# Plot excitatory, inhibitory, and glia contamination markers
VlnPlot(combined, features = c("Slc17a7","Gad1","Snap25","Mbp","Gja1", "Col5a3", "Chodl", "Meis2"))

# Subset glut and gaba ----------------------------------------------------

## In progress ##
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
                     resolution = 0.2)

# doublet clusters 6 and 16,  exclude then recluster
gaba <- subset(gaba, idents = c(as.character(0:5),
                                as.character(7:15),
                                as.character(17:18)))
# recluster 0.2 resolution




VlnPlot(gaba, features = c("Gad1","Kcnc2","Sst","Chodl","Ppp1r1b", "Meis2","Vip","Lamp5", "Ndnf","Cck"))

gaba_ids <- c("Pvalb", "Sst", "Ppp1r1b", "Vip", "Meis2", "Lamp5", "Sst Chodl")
names(gaba_ids) <- levels(gaba)
gaba <- RenameIdents(gaba, gaba_ids)
saveRDS(gaba, file = "gaba.RData")
