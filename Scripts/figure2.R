# Figure 2 plots

# This script generates:
# Figure 2a, excitatory/inhibitory DimPlots
# Figure 2b, excitatory/inhibitory DimPlots grouped by group
# Figure 2c, nuclei prop by group table*
# Figure 2d, DEG table*
# Figure 2e, log2FC for example genes*
# Figure 2f, DEG upset table*

# * denotes that the final plot was made in Prism with output generated in R

# Load packages -----------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(purrr)
library(Libra)
library(rlist)
library(stringr)
library(tibble)
library(DESeq2)
library(UpSetR)
library(ComplexHeatmap)
library(ComplexUpset)
library(data.table)

# Load data ---------------------------------------------------------------

load("glut.RData")
load("gaba.RData")
setwd("~/XPoSE/Output")

# F2a DimPlots ------------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_dimplot.R")

save_dimplot(glut, file_n = "glutclusters.pdf", glutcol = T)

save_dimplot(gaba, file_n = "gabaclusters.pdf", gabacol = T)

# F2b DimPlots and proportions by group ---------------------------------

save_dimplot(glut, groupby = "group", file_n = "glutgroup.pdf", groupcol = T)

save_dimplot(gaba, groupby = "group", file_n = "gabagroup.pdf", groupcol = T)

# F2c proportions by group ----------------------------------------------
source("~/XPoSE/Scripts/Functions/calc_prop.R")

glut_group <- calc_prop(seur_obj = glut, fact1 = 'ratID',
                          fact2 = 'cluster_name',
                          fact3 = 'group',
                          file_n = "glut_group.csv")

gaba_group <- calc_prop(seur_obj = gaba, fact1 = 'ratID',
                          fact2 = 'cluster_name',
                          fact3 = 'group',
                          file_n = "gaba_group.csv")

# F2d DEG table -----------------------------------------------------------

source("~/XPoSE/Scripts/Functions/make_coexp.R")

# Non-active vs Homecage (between subject)
coexp_NH <- make_coexp(seur_obj1 = glut, seur_obj2 = gaba, 
                       threshold = 75, factor = 'group', 
                       comp_vect = c("Non-active","Homecage"))

write.csv(coexp_NH, file = "~/Output/coexp_NH.csv")

# Active vs Homecage (between subject)
coexp_AH <- make_coexp(seur_obj1 = glut, seur_obj2 = gaba, 
                       threshold = 75, factor = 'group', 
                       comp_vect = c("Active","Homecage"))

write.csv(coexp_AH, file = "~/Output/coexp_AH.csv")

# Active vs. Non-active (within subject)
coexp_AN <- make_coexp(seur_obj1 = glut, seur_obj2 = gaba, 
                    threshold = 75, factor = 'group', 
                    comp_vect = c("Active","Non-active"))

write.csv(coexp_AN, file = "~/Output/coexp_AN.csv")

# F2e examples ------------------------------------------------------------

source("~/XPoSE/Scripts/Functions/find_DEcounts.R")

highlights <- c("Vgf", "Ptprn", "Homer1", "Lingo1", "Fosb", "Nptx2", "Scg2",
                "Bdnf","Arc","Pcdh15","Lingo2")

find_DEcounts(ITL56_dds, coexp_AN, "ITL56", highlights)
find_DEcounts(ITL23_dds, coexp_AN, "ITL23", highlights)
find_DEcounts(CTL6_dds, coexp_AN, "CTL6", highlights)
find_DEcounts(PTL5_dds, coexp_AN, "PTL5", highlights)
find_DEcounts(Pvalb_dds, coexp_AN, "Pvalb", highlights)
find_DEcounts(Sst_dds, coexp_AN, "Sst", highlights)

# F2f Upset table ---------------------------------------------------------

source("~/XPoSE/Scripts/Functions/prep_upset.R")
source("~/XPoSE/Scripts/Functions/make_upset.R")

make_upset(prep_upset(coexp_AN), 
           order_vect = c("PTL5","ITL23","ITL56","CTL6","Pvalb","Sst"))