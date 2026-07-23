# Figure 4
library(Seurat)
library(tidyverse)

all <- load('all_annotated_111325.RData')

source('/Users/leej77/Documents/R Work/XPoSE/XPoSE/Scripts/Functions/save_dimplot.R')
source('/Users/leej77/Documents/R Work/XPoSE/XPoSE/Scripts/Functions/calc_prop.R')

# all$experience <- ifelse(all$group == "Homecage", "HC",
#                          "NC")

# Custom groupby-hex combinations
hex_list <- list(
  "cluster_name" = c('CTL6' = '#2D8CB8',
                     'CTL6b' = '#7044AA',
                     'ETL5' = '#0D5A8B',
                     'ITL23' = '#2EBF5E',
                     'ITL5' = '#50B2AD',  
                     'ITL6' = '#58D2CF',
                     'ITvm' = '#B1DE7D',
                     'NPL5' = '#3E9E64',
                     'Pvalb' = '#B9342C',
                     'PvalbChand' = '#FF2D4E',
                     'Sst' = '#FF9900',
                     'SstChodl' = '#B1B10C',
                     'Sncg' = '#D3408D',
                     'Vip' = '#B864CC',
                     'Lamp5' = '#DA808C'),
  "orig.ident" = c('vmPFC1' = '#5A9BC7', 
                   'vmPFC2' = '#89D1D9',
                   'vmPFC3' = '#5A9BC7', 
                   'vmPFC4' = '#89D1D9',
                   'dmPFC1' = '#E08A2D',
                   'dmPFC2' = '#EDCC85',
                   'dmPFC3' = '#E08A2D',
                   'dmPFC4' = '#EDCC85'),
  'sex' = c('male' = "#2C5F2D", 
            'female' = "#97BC62"),
  'experience' =  c('NT' = '#C0C0C0',
                    'RT' = '#AE1E5B',
                    'NC' = '#74b4af',
                    'N' = '#000000'),
  'group' = c('non-active' = '#e37a9e',
              'active' = '#801743')
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
             groupby = "cluster_name",
             splitby = 'experience',
             file_n = 'F4F_all',
             hex_list = hex_list)


# F4G ---------------------------------------------------------------------

seeking <- subset(all, subset = experience == 'RT')

save_dimplot(seeking, 
             groupby = 'group',
             file_n = 'F4G_seeking',
             hex_list = hex_list)

NC <- subset(all, subset = experience == 'NC')

save_dimplot(NC, 
             groupby = 'group',
             file_n = 'F4G_NC',
             hex_list = hex_list)

# F4H ---------------------------------------------------------------------

save_dimplot(seeking, 
             groupby = 'group',
             splitby = 'ratID',
             file_n = 'F4H_seeking',
             hex_list = hex_list)

save_dimplot(NC, 
             groupby = 'group',
             splitby = 'ratID',
             file_n = 'F4H_NC',
             hex_list = hex_list)


# F4I ---------------------------------------------------------------------

# for pie chart
obj_celltype <- calc_prop(seur_obj = all, 
                          fact1 = 'experience',
                          fact2 = 'cluster_name')

write.csv(obj_celltype, file = "F4I_all_clust_by_group.csv")

# for bar plot
clust_prop_cart <- calc_prop(seur_obj = all, 
                             fact1 = 'ratID',
                             fact2 = 'cluster_name',
                             fact3 = 'group')

write.csv(clust_prop_cart, file = "F4I_all_clust_prop_group.csv")


