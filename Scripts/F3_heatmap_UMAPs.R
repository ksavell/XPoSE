library(Seurat)
library(pheatmap)
library(ggplot2)

load("/Users/holmesar/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/Robjects/hc_combined_09112024.RData")

markers <- FindAllMarkers(combine_f)


marker_genes <- c("Slc7a7","Gad1", # general markers
                  "Rfx3", "Cux2", # "ITL23"
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
avg_expr <- AverageExpression(combined_f, features = marker_genes, group.by = "cluster_name")$RNA

expr_matrix <- as.matrix(avg_expr)

expr_matrix <- expr_matrix[marker_genes, cluster_order, drop = FALSE]

# Creating table to show z-scored expression values
scaled_expr_matrix <- t(apply(expr_matrix, 1, scale))
colnames(scaled_expr_matrix) <- colnames(expr_matrix)

write.csv(scaled_expr_matrix,"scaled_expr_matrix.csv", row.names = TRUE) 

# Plot heatmap
pheatmap(scaled_expr_matrix, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#444444", "white", "#000080"))(300),
         main = "Marker Gene Expression by Cluster")

gene_cluster_heatmap <- pheatmap(scaled_expr_matrix, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("#444444", "white", "#000080"))(300),
                                 main = "Marker Gene Expression by Cluster")


# Creating dimplots for excitatory and inhibitory gene expression
          
#legend_colors <- c("#d3d3d3", "#6F77B0", "#06188C")

# Create the DimPlot for the Slc17a
FeaturePlot(combined_f, features = c("Slc17a7"), min.cutoff = 'q1', pt.size = 0.5) +
    theme_void() +   # Removes the background grid
    theme(axis.title = element_blank(),  # Removes axis titles
          legend.position = "right",      # Positions the legend to the right
          #legend.title = element_blank(), # Optional: Remove legend title
          plot.title = element_blank())

Slc17a7_umap <- FeaturePlot(combined_f, features = c("Slc17a7"), min.cutoff = 'q1',pt.size = 0.5) +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "none",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

ggsave("/Users/holmesar/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F3/Slc17a7_fixed.pdf", Slc17a7_umap, width = 10, height = 10)


# Create the DimPlot for the Gad1
FeaturePlot(combined_f, features = c("Gad1"), min.cutoff = 'q1',pt.size = 0.5) +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "right",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

Gad1_umap <- FeaturePlot(combined_f, features = c("Gad1"),min.cutoff = 'q1',pt.size = 0.5) +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = "none",      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())
        
ggsave("/Users/holmesar/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F3/Gad1_fixed.pdf", Gad1_umap, width = 10, height = 10)


# UMAP and gene expression csv generation 
umap_coords <- as.data.frame(combined_f@reductions$umap@cell.embeddings)

slc17a7_expr <- as.data.frame(combined_f@assays$RNA@data["Slc17a7", ])

# Combine UMAP coordinates and gene expression into one data frame
umap_expr_data <- cbind(umap_coords, Slc17a7_Expression = slc17a7_expr)

# Preview the combined data frame
head(umap_expr_data)

write.csv(umap_expr_data, "Slc17a7_UMAP_Data.csv", row.names = TRUE)
