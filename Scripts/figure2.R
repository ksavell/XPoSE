# Figure 2 plots

# This script generates:
# Figure 2a, excitatory/inhibitory DimPlots
# Figure 2b, nuclei counts by group table*
# Figure 2c, excitatory/inhibitory DimPlots grouped by group
# Figure 2d, DEG table*
# Figure 2e, DEG upset table*

# * denotes that the final plot was made in Prism with output generated in R

# Load packages -----------------------------------------------------------
library(Seurat)
library(ggplot2)
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

load("~glut.RData")
load("~gaba.RData")

# F2a DimPlots ------------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_dimplot.R")

save_dimplot(glut, file_n = "glutclusters.pdf", glutcol = T)

save_dimplot(gaba, file_n = "gabaclusters.pdf", gabacol = T)

# F2b cluster percentages by group ----------------------------------------

source("~/XPoSE/Scripts/Functions/calc_counts.R")

calc_counts(seur_obj = glut, fact1 = 'ratID',
            fact2 = 'cluster_name',
            fact3 = 'group', file_n = "group_glut.csv")

calc_counts(seur_obj = gaba, fact1 = 'ratID',
            fact2 = 'cluster_name',
            fact3 = 'group', file_n = "group_gaba.csv")

# F2c DimPlots by group ---------------------------------------------------

save_dimplot(glut, groupby = "group", file_n = "glutgroup.pdf", groupcol = T)

save_dimplot(gaba, groupby = "group", file_n = "gabagroup.pdf", groupcol = T)


# F2d DEG table -----------------------------------------------------------

source("~/XPoSE/Scripts/Functions/make_coexp.R")

coexp_NH <- make_coexp(seur_obj1 = glut, seur_obj2 = gaba, 
                       threshold = 100, factor = 'group', 
                       comp_vect = c("Non-active","Homecage"))

write.csv(coexp_NH, file = "~/Output/coexp_NH.csv")

coexp_AH <- make_coexp(seur_obj1 = glut, seur_obj2 = gaba, 
                       threshold = 100, factor = 'group', 
                       comp_vect = c("Active","Homecage"))

write.csv(coexp_AH, file = "~/Output/coexp_AH.csv")

coexp_AN <- make_coexp(seur_obj1 = glut, seur_obj2 = gaba, 
                    threshold = 100, factor = 'group', 
                    comp_vect = c("Active","Non-active"))

write.csv(coexp_AN, file = "~/Output/coexp_AN.csv")

# F2e Upset table ---------------------------------------------------------

source("~/XPoSE/Scripts/Functions/prep_upset.R")
source("~/XPoSE/Scripts/Functions/make_upset.R")

make_upset(prep_upset(coexp_AN), 
           order_vect = c("IT-L2/3","IT-L5/6","CT-L6","PT-L5","Pvalb","Sst"))

