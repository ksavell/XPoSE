# Figure S2 plots

# This script generates:
# Figure S3a Excitatory Vlnplot
# Figure S2b Inhibitory Vlnplot 
# Figure S2

# * denotes that the final plot was made in Prism with output generated in R

# Load data ---------------------------------------------------------------

load("~glut.RData")
load("~gaba.RData")

# Load packages -----------------------------------------------------------

library(seurat)

# FS2a Exc. Cart QC -------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_vlnplot.R")

save_vlnplot(glut, file_n = "QCglut.pdf", glutcol = T, feature = "nCount_RNA")
save_vlnplot(glut, file_n = "nFeature_glut.pdf", glutcol = T, feature = "nFeature_RNA")

save_vlnplot(gaba, file_n = "QCgaba.pdf", gabacol = T, feature = "nCount_RNA")
save_vlnplot(gaba, file_n = "nFeature_gaba.pdf", gabacol = T, feature = "nFeature_RNA")

