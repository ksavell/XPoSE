# Load necessary packages
library(tidyverse)
library(venn)

get_sig_gene_lists <- function(base_dir, pattern) {
  # Build file pattern and safe filename tag
  file_pattern <- paste0("_group_", pattern, ".csv$")
  output_tag <- gsub("[^A-Za-z0-9]", "_", pattern)
  
  # File matching
  files <- list.files(path = base_dir, pattern = file_pattern, full.names = TRUE)
  
  # Initialize gene lists
  up_gene_lists <- list()
  down_gene_lists <- list()
  
  for (file in files) {
    cluster <- gsub(file_pattern, "", basename(file))
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    
    if (!all(c("gene", "padj", "log2FoldChange") %in% colnames(data))) {
      cat("Skipping", file, "- missing required columns.\n")
      next
    }
    
    up_genes <- unique(na.omit(as.character(data$gene[data$padj < 0.05 & data$log2FoldChange > 0])))
    down_genes <- unique(na.omit(as.character(data$gene[data$padj < 0.05 & data$log2FoldChange < 0])))
    
    if (length(up_genes) > 0) up_gene_lists[[cluster]] <- up_genes
    if (length(down_genes) > 0) down_gene_lists[[cluster]] <- down_genes
    
    if (length(up_genes) == 0 && length(down_genes) == 0) {
      cat("No significant up/downregulated genes in", file, "\n")
    }
  }
  
  # Matrix builder and writer
  write_gene_matrix <- function(gene_list, label) {
    if (length(gene_list) < 2) {
      cat("Not enough clusters with significant genes to create", label, "\n")
      return(NULL)
    }
    
    all_genes <- unique(unlist(gene_list))
    gene_matrix <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
    
    for (cluster in names(gene_list)) {
      gene_matrix[[cluster]] <- ifelse(gene_matrix$Gene %in% gene_list[[cluster]], 1, 0)
    }
    
    out_file <- file.path(base_dir, paste0("group_", output_tag, "_", label, "_overlaps.csv"))
    write.csv(gene_matrix, file = out_file, row.names = FALSE)
    return(gene_matrix)
  }
  
  # Write and return matrices
  up_matrix <- write_gene_matrix(up_gene_lists, "upregulated")
  down_matrix <- write_gene_matrix(down_gene_lists, "downregulated")
  
  # Optional: Save venn diagrams
  make_venn_diagrams <- function(gene_matrix, label) {
    if (is.null(gene_matrix)) return()
    
    # Define clusters of interest
    glut_cols <- c("ITL23", "ITL5", "ITL6", "CTL6", "PTL5")
    gaba_cols <- c("Pvalb", "Sst")
    
    glut_sets <- lapply(glut_cols, function(col) {
      if (col %in% colnames(gene_matrix)) gene_matrix$Gene[gene_matrix[[col]] == 1] else character(0)
    })
    names(glut_sets) <- glut_cols
    
    gaba_sets <- lapply(gaba_cols, function(col) {
      if (col %in% colnames(gene_matrix)) gene_matrix$Gene[gene_matrix[[col]] == 1] else character(0)
    })
    names(gaba_sets) <- gaba_cols
    
    # Save glutamatergic Venn
    pdf(file = file.path(base_dir, paste0(label, "_glut_venn.pdf")))
    venn(glut_sets, ilabels = "counts", zcolor = c(
      'ITL23' = '#41B75F', 
      'ITL5'  = '#5DBFC1', 
      'ITL6'  = '#3A8F87', 
      'CTL6'  = '#2C8CB9', 
      'PTL5'  = '#0A5B8C'
    ))
    dev.off()
    
    # Save GABAergic Venn
    pdf(file = file.path(base_dir, paste0(label, "_gaba_venn.pdf")))
    venn(gaba_sets, ilabels = "counts", zcolor = c(
      'Pvalb' = '#E66027', 
      'Sst'   = '#F8991D'
    ))
    dev.off()
  }
  
  make_venn_diagrams(up_matrix, paste0("group_", output_tag, "_upregulated"))
  make_venn_diagrams(down_matrix, paste0("group_", output_tag, "_downregulated"))
}

# Define working directory
#base_dir <- "~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/Group_Active_HC_2371subset/"

# Active v homecage
get_sig_gene_lists(base_dir = "/Users/savellke/Projects/XPoSE/DES_Active_Homecage_equalratgroup", pattern = "Active_Homecage")

# Active v Nonactive
get_sig_gene_lists(base_dir = "/Users/savellke/Projects/XPoSE/DES_Active_Nonactive_equalratgroup", pattern = "Active_Non-active")

