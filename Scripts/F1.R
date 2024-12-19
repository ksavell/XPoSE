# Figure 1

# Info --------------------------------------------------------------------

# This script 


# Load libraries and data -------------------------------------------------

library(Seurat)
library(tidyverse)

# Count neuron/non-neuron by rat/capture  ---------------------------------

load("combined_withmmmannotation_09112024.RData")

source("Scripts/Functions/calc_prop.R")

obj_celltype <- calc_prop(seur_obj = combined, 
                          fact1 = 'ratID',
                          fact2 = 'celltype',
                          fact3 = 'orig.ident')

write.csv(obj_celltype, file = "celltype_by_rat_cart.csv")


# HC exp DimPlots ---------------------------------------------------------

load("hc_combined_09112024.RData")

hc <- combined_f

source("Scripts/Functions/save_dimplot.R")

# Custom groupby-hex combinations
hex_list <- list(
  'cluster_name' = c('CTL6' = '#2C8CB9',
                     'PTL5' = '#0A5B8C',
                     'ITL23' = '#41B75F',
                     'ITL5' = '#5DBFC1',  
                     'ITL6' = '#3A8F87',
                     'NPL56' = '#3C9E64',
                     'CTL6b' = '#6F499D',
                     'Pvalb' = '#E66027',
                     'Sst' = '#F8991D',
                     'Meis2' = '#C52126',
                     'Vip' = '#A669AB',
                     'Lamp5' = '#DB808C',
                     'SstChodl' = '#B0B235',
                     'PvalbChand' = '#AD6C49'),
  'orig.ident' = c('C1' = '#22677A', 
                   'C2' = '#D4A841'),
  'sex' = c('male' = "#2C5F2D", 
            'female' = "#97BC62")
)

save_dimplot(hc, 
             groupby = 'cluster_name',
             file_n = 'hc',
             hex_list = hex_list)

save_dimplot(hc, 
             groupby = 'orig.ident',
             splitby = 'ratID',
             file_n = 'hc',
             hex_list = hex_list)


# Cluster prop by rat/capture ---------------------------------------------

clust_prop_cart <- calc_prop(seur_obj = hc, 
                          fact1 = 'ratID',
                          fact2 = 'cluster_name',
                          fact3 = 'capture') 

write.csv(clust_prop_cart, file = "clust_prop_cart.csv")


