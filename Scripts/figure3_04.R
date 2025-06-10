# Figure 3

# This script generates data for all panels in F3

# Loading -----------------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)

source("Scripts/Functions/save_dimplot.R")
source("Scripts/Functions/calc_prop.R")

# load in clustered hc object that is output of createobject_01.R

load("hc_09112024.RData")
hc <- combined_f

# define colors
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

# F3A ---------------------------------------------------------------------

save_dimplot(hc, 
             groupby = 'cluster_name',
             file_n = 'hc',
             hex_list = hex_list)


# F3B ---------------------------------------------------------------------

# feature plot save

# Create the DimPlot for the Slc17a7
FeaturePlot(hc, features = c("Slc17a7"), 
            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

Slc17a7_umap <- FeaturePlot(hc, features = c("Slc17a7"), 
                            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

ggsave("/Users/holmesar/Desktop/Slc17a7.svg", Slc17a7_umap, width = 10, height = 10)


# Create the DimPlot for the Gad1
FeaturePlot(hc, features = c("Gad1"), 
            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

Gad1_umap <- FeaturePlot(hc, features = c("Gad1"), 
                         cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

ggsave("/Users/holmesar/Desktop/Gad1.svg", Gad1_umap, width = 10, height = 10)

# F3C ---------------------------------------------------------------------

# marker heatmap

marker_genes <- c("Rfx3", "Cux2", # "ITL23"
                  "Rorb", "Slc7a11", # "ITL5"
                  "Col6a1", "Col6a2", # "ITL6"
                  "Syt6", "Foxp2", # "CTL6"
                  "Ctgf", "Cplx3", # "CTL6b"
                  "Gpc5", "Fezf2", # "PTL5"
                  "Tshz2", "Htr4", # "NPL56"
                  "Erbb4", "Kcnc2", # "Pvalb"
                  "Sst", "Elfn1", # "Sst"
                  "Slc6a1", "Unc5b", # "PvalbChand"
                  "Chodl", "Nos1", # "SstChodl"
                  "Vip", "Adarb2", # "Vip"
                  "Lamp5", "Egfr", # "Lamp5"
                  "Meis2", "Drd1", "Drd2")  # "Meis2"

cluster_order <- c("ITL23", "ITL5", "ITL6", "CTL6", "CTL6b", "PTL5", "NPL56", "Pvalb", "Sst", "PvalbChand", "SstChodl", "Vip", "Lamp5", "Meis2")

# Compute average expression per cluster
avg_expr <- AverageExpression(hc, features = marker_genes, group.by = "cluster_name")$RNA

expr_matrix <- as.matrix(avg_expr)

expr_matrix <- expr_matrix[marker_genes, cluster_order, drop = FALSE]

# Creating table to show z-scored expression values
scaled_expr_matrix <- t(apply(expr_matrix, 1, scale))
colnames(scaled_expr_matrix) <- colnames(expr_matrix)

write.csv(scaled_expr_matrix,"scaled_expr_matrix.csv", row.names = TRUE) 


# F3D ---------------------------------------------------------------------

save_dimplot(hc, 
             groupby = 'orig.ident',
             splitby = 'ratID',
             file_n = 'hc',
             hex_list = hex_list)


# F3E ---------------------------------------------------------------------

# props 

clust_prop_cart <- calc_prop(seur_obj = hc, 
                             fact1 = 'ratID',
                             fact2 = 'cluster_name',
                             fact3 = 'orig.ident') 

write.csv(clust_prop_cart, file = "f3e_clust_prop_cart.csv")