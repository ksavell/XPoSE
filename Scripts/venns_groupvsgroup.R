# venns for AvsHC and AvsNA

# Load required library
library(venn)

# Define working directory
base_dir <- "~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs"

# Get lists of files for Active_Homecage and Active_Nonactive
homecage_files <- list.files(path = base_dir, pattern = "_group_Active_Homecage.csv$", full.names = TRUE)
nonactive_files <- list.files(path = base_dir, pattern = "_group_Active_Nonactive.csv$", full.names = TRUE)

# Extract cluster names
clusters <- gsub("_group_Active_Homecage.csv", "", basename(homecage_files))

# Initialize lists to store significant genes for each condition
homecage_genes <- list()
nonactive_genes <- list()

# Process each cluster file
for (cluster in clusters) {
  homecage_file <- file.path(base_dir, paste0(cluster, "_group_Active_Homecage.csv"))
  nonactive_file <- file.path(base_dir, paste0(cluster, "_group_Active_Nonactive.csv"))
  
  if (file.exists(homecage_file)) {
    data_homecage <- read.csv(homecage_file, header = TRUE, stringsAsFactors = FALSE)
    if (all(c("gene", "padj") %in% colnames(data_homecage))) {
      homecage_genes[[cluster]] <- unique(na.omit(data_homecage$gene[data_homecage$padj < 0.05]))
    }
  }
  
  if (file.exists(nonactive_file)) {
    data_nonactive <- read.csv(nonactive_file, header = TRUE, stringsAsFactors = FALSE)
    if (all(c("gene", "padj") %in% colnames(data_nonactive))) {
      nonactive_genes[[cluster]] <- unique(na.omit(data_nonactive$gene[data_nonactive$padj < 0.05]))
    }
  }
}

# Create Venn diagrams for each cluster using the venn package
pdf(file = file.path(base_dir, "Active_vs_homecage_or_nonactive_venn.pdf"))

for (cluster in clusters) {
  if (length(homecage_genes[[cluster]]) > 0 || length(nonactive_genes[[cluster]]) > 0) {
    venn_data <- list(
      "Active Homecage" = homecage_genes[[cluster]],
      "Active Nonactive" = nonactive_genes[[cluster]]
    )
    
    venn(venn_data, ilabels = "counts", zcolor = c("blue", "red"), main = paste("Cluster:", cluster))
  }
}

dev.off()
