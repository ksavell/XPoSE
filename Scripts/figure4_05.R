# Figure 4
library(Seurat)
library(tidyverse)


source('/Users/leej77/Documents/R Work/XPoSE/XPoSE/Scripts/Functions/save_dimplot.R')
source('/Users/leej77/Documents/R Work/XPoSE/XPoSE/Scripts/Functions/calc_prop.R')

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
  "orig.ident" = c('C1' = '#5A9BC7',
                   'C2' = '#E08A2D'),
  'sex' = c('male' = "#2C5F2D", 
            'female' = "#97BC62"),
  'experience' =  c('HC' = '#C0C0C0',
                    'RT' = '#AE1E5B',
                    'NC' = '#0A8C81',
                    'N' = 'black'),
  'group' = c('non-active' = '#75C3BC',
              'active' = '#00675F')
)

# F4D ---------------------------------------------------------------------

# post QC NC counts

obj_celltype <- calc_prop(seur_obj = obj, 
                          fact1 = 'ratID',
                          fact2 = 'group')

write.csv(obj_celltype, file = "f4d_group_by_rat.csv")

# F4E ---------------------------------------------------------------------

# Dimplots for HC/NC combined object

save_dimplot(obj, 
             groupby = 'cluster_name',
             file_n = 'F4E_all',
             hex_list = hex_list)


# F4F ---------------------------------------------------------------------

save_dimplot(obj, 
             groupby = "experience",
             file_n = 'F4F_all',
             hex_list = hex_list)


# F4G ---------------------------------------------------------------------

NC <- subset(obj, subset = experience == 'NC')

save_dimplot(NC, 
             groupby = 'group',
             file_n = 'F4G_NC',
             hex_list = hex_list)

# F4H ---------------------------------------------------------------------

save_dimplot(NC, 
             groupby = 'group',
             splitby = 'ratID',
             file_n = 'F4H_NC',
             hex_list = hex_list)


# F4I ---------------------------------------------------------------------

# for pie chart
obj_celltype <- calc_prop(seur_obj = obj, 
                          fact1 = 'experience',
                          fact2 = 'cluster_name')

write.csv(obj_celltype, file = "F4I_all_clust_by_group.csv")

# for bar plot
clust_prop_cart <- calc_prop(seur_obj = obj, 
                             fact1 = 'ratID',
                             fact2 = 'cluster_name',
                             fact3 = 'group')

write.csv(clust_prop_cart, file = "F4I_all_clust_prop_group.csv")

all_p_values <- c(0.138772,
                  0.095136,
                  0.882892,
                  0.619633,
                  0.989023,
                  0.048194,
                  0.000008,
                  0.457147, # this ends NC:NA v HC
                  0.001121,
                  0.165038,
                  0.000180,
                  0.040186,
                  0.000001,
                  0.002378,
                  0.000599,
                  0.000159, # this ends NC:A vs HC
                  0.000107,
                  0.094900,
                  0.003804,
                  0.028463,
                  0.000214,
                  0.001804,
                  0.002276,
                  0.010723) # this ends NC:A vs NC:NA


adjusted_p_values <- p.adjust(all_p_values, method = "BH")

write.csv(adjusted_p_values, "fig4i_adjpval_BH.csv")

write.csv(clustgroup_prop_04282026, "cluster_proportions_05052026.csv")

