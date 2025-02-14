# Load necessary packages
library(dplyr)
library(ggvenn)
library(venn)

# Define working directory
base_dir <- "~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/Group_Active_HC_2371subset/"

# Get list of all files in the directory
files <- list.files(path = base_dir, pattern = "_group_Active_Homecage.csv$", full.names = TRUE)

# Initialize a list to store significant genes for each cluster
sig_gene_lists <- list()

# Process each file
for (file in files) {
  cluster <- gsub("_group_Active_Homecage.csv", "", basename(file))  # Extract cluster name
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  if (!all(c("gene", "padj") %in% colnames(data))) {
    cat("Skipping", file, "- missing required columns.\n")
    next
  }
  
  sig_genes <- unique(as.character(na.omit(data$gene[data$padj < 0.05])))
  
  if (length(sig_genes) > 0) {
    sig_gene_lists[[cluster]] <- sig_genes
  } else {
    cat("No significant genes in", file, "\n")
  }
}

# Ensure at least two clusters have significant genes
if (length(sig_gene_lists) < 2) {
  cat("Not enough clusters with significant genes to create an UpSet plot.\n")
} else {
  # Create a presence/absence matrix for all genes across clusters
  all_genes <- unique(unlist(sig_gene_lists))
  gene_matrix <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
  
  for (cluster in names(sig_gene_lists)) {
    gene_matrix[[cluster]] <- ifelse(gene_matrix$Gene %in% sig_gene_lists[[cluster]], 1, 0)
  }
  
  # Save the matrix as a CSV
  write.csv(gene_matrix, file = file.path(base_dir, "group_Active_Homecage_overlaps.csv"), row.names = FALSE)
}

# --- Convert Binary Matrix to List of Genes for `venn()` ---
# Define the subsets of columns for "glut" and "gaba"
glut_cols <- c("ITL23","ITL5","ITL6","CTL6", "PTL5")  #  glutamatergic clusters
gaba_cols <- c("Pvalb","Sst")  #  GABAergic clusters

# Convert for Glut Venn
glut_sets <- lapply(glut_cols, function(col) gene_matrix$Gene[gene_matrix[[col]] == 1])
names(glut_sets) <- glut_cols

# Convert for GABA Venn
gaba_sets <- lapply(gaba_cols, function(col) gene_matrix$Gene[gene_matrix[[col]] == 1])
names(gaba_sets) <- gaba_cols

# Create new dataframe with 'Gene' column and collapsed 'glut' and 'gaba' categories
glut_gaba_sets <- gene_matrix %>%
  mutate(
    glut = ifelse(rowSums(select(., all_of(glut_cols))) > 0, 1, 0),
    gaba = ifelse(rowSums(select(., all_of(gaba_cols))) > 0, 1, 0)
  ) %>%
  select(Gene, glut, gaba)  # Keep only relevant columns

write.csv(glut_gaba_sets, "glut_gaba_geneoverlaps.csv")
glut_gaba_sets <- glut_gaba_sets %>%
  select(glut,gaba)

# --- Save Venn Diagrams as PDFs ---
pdf(file = file.path(base_dir, "glut_venn.pdf"))
venn(glut_sets, ilabels = "counts", zcolor = c('ITL23' = '#41B75F', 
                                               'ITL5' = '#5DBFC1', 
                                               'ITL6' = '#3A8F87', 
                                               'CTL6' = '#2C8CB9', 
                                               'PTL5' = '#0A5B8C'))
dev.off()

pdf(file = file.path(base_dir, "gaba_venn.pdf"))
venn(gaba_sets, ilabels = "counts", zcolor = c('Pvalb' = '#E66027', 'Sst' = '#F8991D'))
dev.off()

pdf(file = file.path(base_dir, "glut_vs_gaba_venn.pdf"))
venn(glut_gaba_sets, ilabels = "counts",
     zcolor = c('#1E9BB5','#D35400'))
dev.off()


