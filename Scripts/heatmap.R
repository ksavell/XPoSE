library(Seurat)
library(pheatmap)

# Define your marker genes
marker_genes <- c("Slc17a7", "Gad1", "Sst", "Rbx3")  # Replace with your actual genes

# Compute average expression per cluster
avg_expr <- AverageExpression(all, features = marker_genes, group.by = "cluster_name")$RNA

# Convert to a matrix for visualization
expr_matrix <- as.matrix(avg_expr)

# Plot heatmap
pheatmap(expr_matrix, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Marker Gene Expression by Cluster")
