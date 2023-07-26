# Figure 2 plots

# This script generates:
# Figure 2a, excitatory/inhibitory DimPlots
# Figure 2b, cluster percentages by group table*
# Figure 2c, excitatory/inhibitory DimPlots grouped by group
# Figure 2d, DEG table*
# Figure 2e, DEG upset table*

# * denotes that the final plot was made in Prism from output generated in R

# Load packages -----------------------------------------------------------
library(Seurat)
library(ggplot2)

# Load data ---------------------------------------------------------------

# set working directory
setwd("~/Documents/GitHub/XPoSE/Scripts/Output/")

load("~/Input/glut.RData")
load("~/Input/gaba.RData")


# F2a DimPlots ------------------------------------------------------------

source("~/XPoSE/Scripts/Functions/save_dimplot.R")

save_dimplot(glut, savename = "glutclusters") ##this isn't working yet

