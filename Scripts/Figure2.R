# Figure 2 plots

# This script generates:
# Figure 2a, excitatory/inhibitory DimPlots
# Figure 2b, cluster percentages by group table*
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

# in progress

# F2c DimPlots by group ---------------------------------------------------

save_dimplot(glut, groupby = "group", file_n = "glutgroup.pdf", groupcol = T)

save_dimplot(gaba, groupby = "group", file_n = "gabagroup.pdf", groupcol = T)


# F2d DEG table -----------------------------------------------------------

coexp_PN <- make_coexp(glut, gaba, 100,  "group", c("positive", "negative"))
coexp_NH <- make_coexp(glut, gaba, 100,  "group", c("negative", "homecage"))
coexp_PN <- make_coexp(glut, gaba, 100,  "group", c("positive", "homecage"))


# F2e Upset table ---------------------------------------------------------

make_upset(prep_upset(coexp_PN), c("CT L6", "IT L2/3", "IT L5/6", "PT L5", "Pvalb", "Sst"))
