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

# Load data ---------------------------------------------------------------

load("~/Input/glut.RData")
load("~/Input/gaba.RData")


# F2a DimPlots ------------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_dimplot.R")

save_dimplot(glut, file_n = "glutclusters.pdf")

save_dimplot(gaba, file_n = "gabaclusters.pdf")


# F2b cluster percentages by group ----------------------------------------

# in progress


# F2c DimPlots by group ---------------------------------------------------

save_dimplot(glut, groupby = "group", file_n = "glutgroup.pdf")

save_dimplot(gaba, groupby = "group", file_n = "gabagroup.pdf")


# F2d DEG table -----------------------------------------------------------

# in progress, DJT! 


# F2e Upset table ---------------------------------------------------------

# in progress, DJT!
