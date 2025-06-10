# Figure 4

library(Seurat)
library(tidyverse)

load("~/Projects/XPoSE/all_10312024.RData")

source("Scripts/Functions/calc_prop.R")
source("Scripts/Functions/save_dimplot.R")

all$experience <- ifelse(all$group == "Homecage", "HC",
                         "NC")

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
            'female' = "#97BC62"),
  'experience' =  c('HC' = '#c0c0c0',
                    'NC' = '#ae1e5b'),
  'group' = c('Non-active' = '#e37a9e',
              'Active' = '#801743')
)

# F4D ---------------------------------------------------------------------

# post QC NC counts

obj_celltype <- calc_prop(seur_obj = all, 
                          fact1 = 'ratID',
                          fact2 = 'group')

write.csv(obj_celltype, file = "f4d_group_by_rat.csv")

# F4E ---------------------------------------------------------------------

# Dimplots for HC/NC combined object

save_dimplot(all, 
             groupby = 'cluster_name',
             file_n = 'F4E_all',
             hex_list = hex_list)


# F4F ---------------------------------------------------------------------

save_dimplot(all, 
             groupby = 'experience',
             file_n = 'F4F_all',
             hex_list = hex_list)


# F4G ---------------------------------------------------------------------

nc <- subset(all, subset = experience == 'NC')

save_dimplot(nc, 
             groupby = 'group',
             file_n = 'F4G_nc',
             hex_list = hex_list)

# F4H ---------------------------------------------------------------------

save_dimplot(nc, 
             groupby = 'group',
             splitby = 'ratID',
             file_n = 'F4H_nc',
             hex_list = hex_list)


# F4I ---------------------------------------------------------------------

# for pie chart
obj_celltype <- calc_prop(seur_obj = all, 
                          fact1 = 'group',
                          fact2 = 'cluster_name')

write.csv(obj_celltype, file = "F4I_all_clust_by_group.csv")

# for bar plot
clust_prop_cart <- calc_prop(seur_obj = all, 
                             fact1 = 'ratID',
                             fact2 = 'cluster_name',
                             fact3 = 'group')

write.csv(clust_prop_cart, file = "F4I_all_clust_prop_group.csv")