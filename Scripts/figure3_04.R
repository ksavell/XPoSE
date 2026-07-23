# Figure 3
# This script generates data for all panels in F3

# Loading -----------------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)

source('/Users/leej77/Documents/R Work/XPoSE/XPoSE/Scripts/Functions/save_dimplot.R')
source('/Users/leej77/Documents/R Work/XPoSE/XPoSE/Scripts/Functions/calc_prop.R')

# load in clustered hc object that is output of createobject_01.R

all <- load("objects/dmvmPFC_annotated07162026")

# define colors
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
  "sex" = c('male' = "#2C5F2D", 
            'female' = "#97BC62")
)

# F3A ---------------------------------------------------------------------

save_dimplot(all, 
             groupby = 'cluster_name',
             file_n = 'all',
             hex_list = hex_list)


# F3B ---------------------------------------------------------------------

# feature plot save
legend_colors <- c('#D1D1D1', '#2D00FF')
# Create the DimPlot for the Slc17a7
FeaturePlot(all, features = c("Slc17a7"), 
            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

Slc17a7_umap <- FeaturePlot(all, features = c("Slc17a7"), 
                            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

ggsave("Slc17a7.svg", Slc17a7_umap, width = 10, height = 10)


# Create the DimPlot for the Gad1
FeaturePlot(all, features = c("Gad1"), 
            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

Gad1_umap <- FeaturePlot(all, features = c("Gad1"), 
                         cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

ggsave("Gad1.svg", Gad1_umap, width = 10, height = 10)

# F3C ---------------------------------------------------------------------

# marker heatmap

marker_genes <- c("Rfx3", "Cux2", # "ITL23"
                  "Rorb", "Slc7a11", # "ITL5"
                  "Col6a1", "Col6a2", # "ITL6"
                  "Ndst4", "Nrp2",# "ITvm"
                  "Syt6", "Foxp2", # "CTL6"
                  "Ctgf", "Cplx3", # "CTL6b"
                  "Gpc5", "Fezf2", # "ETL5"
                  "Tshz2", "Htr4", # "NPL5"
                  "F2r", "Kcnc2", # "Pvalb"
                  "Sst", "Elfn1", # "Sst"
                  "Slc6a1", "Unc5b", # "PvalbChand"
                  "Chodl", "Nos1", # "SstChodl"
                  "Vip", "Prox1", # "Vip"
                  "Lamp5", "Egfr", # "Lamp5"
                  "Htr3a", "Frem1" # "Sncg"
                  )

cluster_order <- c("ITL23", "ITL5", "ITL6", "ITvm", "CTL6", "CTL6b", "ETL5", "NPL5", "Pvalb", "Sst", "PvalbChand", "SstChodl", "Vip", "Lamp5", "Sncg")

# Compute average expression per cluster
marker_genes <- marker_genes
avg_expr <- AverageExpression(all, features = marker_genes, group.by = "cluster_name")$RNA

expr_matrix <- as.matrix(avg_expr)

expr_matrix <- expr_matrix[marker_genes, cluster_order, drop = FALSE]

# Creating table to show z-scored expression values
scaled_expr_matrix <- t(apply(expr_matrix, 1, scale))
colnames(scaled_expr_matrix) <- colnames(expr_matrix)

write.csv(scaled_expr_matrix,"scaled_expr_matrix.csv", row.names = TRUE) 


# F3D ---------------------------------------------------------------------

save_dimplot(all, 
             groupby = "orig.ident",
             splitby = "ratID",
             file_n = "all",
             hex_list = hex_list)

# F3E ---------------------------------------------------------------------

# props 

clust_prop_cart <- calc_prop(seur_obj = all, 
                             fact1 = 'ratID',
                             fact2 = 'cluster_name',
                             fact3 = 'orig.ident') 

write.csv(clust_prop_cart, file = "f3e_clust_prop_cart.csv")

