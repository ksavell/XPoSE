# Load necessary libraries
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Cairo)

# Define the working directory (update if needed)
base_dir <- "~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs"

# Get list of all files in the directory
files <- list.files(path = base_dir, pattern = "_group_Active_Homecage.csv$", full.names = TRUE)

# Initialize a list to store significant genes for each cluster
sig_gene_lists <- list()

# Process each file
for (file in files) {
  # Extract cluster name from the filename
  cluster <- gsub("_group_Active_Homecage.csv", "", basename(file))
  
  # Read the file
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Check if required columns exist
  if (!all(c("gene", "padj") %in% colnames(data))) {
    cat("Skipping", file, "- missing required columns.\n")
    next
  }
  
  # Remove rows with missing values in gene or padj
  data <- na.omit(data[, c("gene", "padj")])
  
  # Filter for significant genes (padj < 0.05)
  sig_genes <- unique(as.character(data$gene[data$padj < 0.05]))
  
  # Store only if there are significant genes
  if (length(sig_genes) > 0) {
    sig_gene_lists[[cluster]] <- sig_genes
  } else {
    cat("No significant genes in", file, "\n")
  }
}

# Ensure we have at least two clusters with significant genes
if (length(sig_gene_lists) < 2) {
  cat("Not enough clusters with significant genes to create an UpSet plot.\n")
} else 
  {
  # Create a presence/absence matrix for all genes across clusters
  all_genes <- unique(unlist(sig_gene_lists))  # Get all unique genes
  gene_matrix <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
  
  # Fill in presence/absence for each cluster
  for (cluster in names(sig_gene_lists)) {
    gene_matrix[[cluster]] <- ifelse(gene_matrix$Gene %in% sig_gene_lists[[cluster]], 1, 0)
  }
  
  # Remove the gene column to get a binary matrix
  gene_matrix_binary <- gene_matrix %>% select(-Gene)
  
  # Debugging: Check if the matrix is correctly populated
  print(dim(gene_matrix_binary))  # Check dimensions
  print(head(gene_matrix_binary))  # Print first few rows
  
  # Ensure gene_matrix has at least one non-zero entry before plotting
  if (all(gene_matrix_binary == 0)) {
    stop("Error: All values in the gene_matrix are zero! No overlaps to plot.")
  }
  
  # Create an UpSet plot
  # CairoPDF("UpSet_Active_Homecage.pdf", width = 20, height = 10)
  # upset(
  #   gene_matrix_binary, 
  #   sets = names(sig_gene_lists), 
  #   keep.order = TRUE,
  #   sets.bar.color = "steelblue", 
  #   order.by = "freq", 
  #   mainbar.y.label = "Number of Shared Significant Genes",
  #   sets.x.label = "Number of Genes per Cluster"
  # )
  # dev.off()
  
  p <- upset(
    gene_matrix_binary, 
    sets = names(cluster), 
    keep.order = TRUE,
    sets.bar.color = "steelblue", 
    order.by = "freq", 
    mainbar.y.label = "Number of Shared Significant Genes",
    sets.x.label = "Number of Genes per Cluster"
  )
  
  ggsave("UpSet_Active_Homecage.pdf", plot = p, width = 20, height = 10)
  
  write.csv(gene_matrix, file = "group_Active_Homecage_overlaps.csv")
}


