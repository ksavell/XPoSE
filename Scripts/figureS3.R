# Figure S3 plots

# This script generates:
# Figure S3a Genes Vlnplot
# Figure S3b Transcripts Vlnplot 
# Figure S2c Heatmaps

# * denotes that the final plot was made in Prism with output generated in R

# Load data ---------------------------------------------------------------

load("~glut.RData")
load("~gaba.RData")

# Load packages -----------------------------------------------------------

library(seurat)
library(ggplot2)

# FS3a Gene VlnPlots ------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_vlnplot.R")

save_vlnplot(glut, file_n = "nFeature_glut.pdf", 
             glutcol = T, feature = "nFeature_RNA")
save_vlnplot(gaba, file_n = "nFeature_gaba.pdf", 
             gabacol = T, feature = "nFeature_RNA")

# FS3b Transcript VlnPlots ------------------------------------------------

save_vlnplot(glut, file_n = "nCount_glut.pdf", 
             glutcol = T, feature = "nCount_RNA")
save_vlnplot(gaba, file_n = "nCount_gaba.pdf", 
             gabacol = T, feature = "nCount_RNA")

# FS3c Heatmaps -----------------------------------------------------------

source("~/XPoSE/Scripts/Functions/make_heatmap.R")

make_heatmap(glut, groupcol = c("#64C7C8",'#41B75F','#2C8CB9','#0A5B8C',
                                '#3C9E64','#6F499D'), 
             pdffilen = "heatmap_glut",
             csvfilen = "glutmarkers")

make_heatmap(gaba, groupcol = c('#E66027','#F8991D',
                                '#C03C82','#A669AB',
                                '#C52126','#DB808C',
                                '#B0B235'), 
             pdffilen = "heatmap_glut",
             csvfilen = "glutmarkers")
