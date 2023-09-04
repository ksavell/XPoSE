# Figure S1 plots

# This script generates:
# Figure S1a Excitatory UMAP split by cartridge
# Figure S1b Inhibitory UMAP split by cartridge

# * denotes that the final plot was made in Prism with output generated in R

# Load data ---------------------------------------------------------------

setwd("~/XPoSE")
load("glut.RData")
load("gaba.RData")

library(Seurat)
library(ggplot2)

# FS1a Exc. Cart QC -------------------------------------------------------

setwd("~/XPoSE/Output")
source("~/XPoSE/Scripts/Functions/save_dimplot.R")

save_dimplot(glut, file_n = "glutcart.pdf", glutcol = T, 
             splitby = "orig.ident", splitl = 2)

# FS1b Inh. Cart QC -------------------------------------------------------

save_dimplot(gaba, file_n = "gabacart.pdf", gabacol = T, 
             splitby = "orig.ident", splitl = 2)
