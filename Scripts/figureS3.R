# Figure S3 plots

# This script generates:
# Figure S1a Excitatory UMAP split by sex, percent plot*
# Figure S1b Inhibitory UMAP split by sex, percent plot*

# * denotes that the final plot was made in Prism with output generated in R

# Load data ---------------------------------------------------------------

load("~glut.RData")
load("~gaba.RData")

# Load packages -----------------------------------------------------------

library(seurat)
library(ggplot2)

# FS3a Exc. by sex --------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_dimplot.R")

save_dimplot(glut, file_n = "glutsec.pdf", glutcol = T, 
             splitby = "sex", splitl = 2)

source("~/XPoSE/Scripts/Functions/calc_counts.R")

calc_counts(seur_obj = glut, fact1 = 'ratID',
            fact2 = 'seurat_clusters',
            fact3 = 'sex', file_n = "sex_glut.csv")

# FS1b Inh. by sex --------------------------------------------------------

save_dimplot(gaba, file_n = "gabasex.pdf", gabacol = T, 
             splitby = "sex", splitl = 2)

calc_counts(seur_obj = gaba, fact1 = 'ratID',
            fact2 = 'seurat_clusters',
            fact3 = 'sex', file_n = "sex_gaba.csv")
