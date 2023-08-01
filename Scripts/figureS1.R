# Figure S1 plots

# This script generates:
# Figure S1a Excitatory UMAP split by cartridge, percent plot*
# Figure S1b Inhibitory UMAP split by cartridge, percent plot*

# * denotes that the final plot was made in Prism with output generated in R

# Load data ---------------------------------------------------------------

load("~glut.RData")
load("~gaba.RData")

library(Seurat)
library(ggplot2)

# FS1b Exc. Cart QC -------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_dimplot.R")

save_dimplot(glut, file_n = "glutcart.pdf", glutcol = T, 
             splitby = "orig.ident", splitl = 2)

source("~/XPoSE/Scripts/Functions/calc_counts.R")

calc_counts(seur_obj = glut, fact1 = 'ratID',
                 fact2 = 'seurat_clusters',
                 fact3 = 'orig.ident', file_n = "cart_glut.csv")

# FS1b Inh. Cart QC -------------------------------------------------------

save_dimplot(gaba, file_n = "gabacart.pdf", gabacol = T, 
             splitby = "orig.ident", splitl = 2)

calc_counts(seur_obj = gaba, fact1 = 'ratID',
            fact2 = 'seurat_clusters',
            fact3 = 'orig.ident', file_n = "cart_gaba.csv")

