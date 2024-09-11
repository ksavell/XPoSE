# Figure 1

# Info --------------------------------------------------------------------

# This script 


# Load libraries and data -------------------------------------------------

library(Seurat)
library(tidyverse)

load("combined09112024.RData")


# Count neuron/non-neuron by rat/capture  ---------------------------------

source("Scripts/Functions/calc_prop.R")

obj_celltype <- calc_prop(seur_obj = combined, 
                          fact1 = 'ratID',
                          fact2 = 'celltype',
                          fact3 = 'orig.ident')

write.csv(obj_celltype, file = "celltype_by_rat_cart.csv")
