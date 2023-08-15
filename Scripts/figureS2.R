# Figure S2 plots

# This script generates:
# Figure S2a Genes Vlnplot
# Figure S2b Transcripts Vlnplot 
# Figure S2c Heatmaps

# * denotes that the final plot was made in Prism with output generated in R

# Load data ---------------------------------------------------------------

load("~glut.RData")
load("~gaba.RData")

# Load packages -----------------------------------------------------------

library(seurat)
library(ggplot2)
library(dplyr)

# FS2a Gene VlnPlots ------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_vlnplot.R")

save_vlnplot(glut, file_n = "nFeature_glut.pdf", 
             glutcol = T, feature = "nFeature_RNA",
             vln_max = 5000)
save_vlnplot(gaba, file_n = "nFeature_gaba.pdf", 
             gabacol = T, feature = "nFeature_RNA",
             vln_max = 5000)

# FS2b Transcript VlnPlots ------------------------------------------------

save_vlnplot(glut, file_n = "nCount_glut.pdf", 
             glutcol = T, feature = "nCount_RNA",
             vln_max = 15000)
save_vlnplot(gaba, file_n = "nCount_gaba.pdf", 
             gabacol = T, feature = "nCount_RNA",
             vln_max = 15000)

# FS2c Heatmaps -----------------------------------------------------------

source("~/XPoSE/Scripts/Functions/make_heatmap.R")

make_heatmap(glut, groupcol = c('#0A5B8C','#64C7C8',
                                '#2C8CB9','#41B75F',
                                '#6F499D','#3C9E64'), 
             pdffilen = "heatmap_glut",
             csvfilen = "glutmarkers")

make_heatmap(gaba, groupcol = c('#F8991D','#E66027',
                                '#C03C82','#C52126',
                                '#B0B235','#DB808C',
                                '#A669AB'), 
             pdffilen = "heatmap_gaba",
             csvfilen = "gabamarkers")
