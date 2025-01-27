# overlap in consistent genes

# Load necessary library
library(VennDiagram)
library(grid)

# Define the directory containing all cluster folders
base_dir <- "Droprat"  # Replace with the path to your base directory
output_file <- "Consistent_unique_shared_counts.csv"  # Output CSV file

# Initialize a data frame to store results
results <- data.frame(
  Cluster = character(),
  Uniquely_DESeq2 = integer(),
  Uniquely_FindMarkers = integer(),
  Shared = integer(),
  stringsAsFactors = FALSE
)

# Get a list of cluster directories
cluster_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)


# Loop through each cluster
for (cluster_dir in cluster_dirs) {
  # Extract cluster ID from the directory name
  cluster_id <- basename(cluster_dir)
  
  # Define file paths for the cluster
  deseq2_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnConsistent_Group_Active_Homecage_", cluster_id, ".csv"))
  findmarkers_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnConsistent_Group_Active_Homecage_", cluster_id, ".csv"))
  
  # Check if both files exist
  if (file.exists(deseq2_file) && file.exists(findmarkers_file)) {
    # Read data
    deseq2_data <- read.csv(deseq2_file, header = TRUE, stringsAsFactors = FALSE)
    findmarkers_data <- read.csv(findmarkers_file, header = TRUE, stringsAsFactors = FALSE)
    
    # Extract gene names (assuming the first column is gene names)
    deseq2_genes <- na.omit(deseq2_data[[1]])
    findmarkers_genes <- na.omit(findmarkers_data[[1]])
    
    # Calculate unique and shared genes
    unique_deseq2 <- setdiff(deseq2_genes, findmarkers_genes)
    unique_findmarkers <- setdiff(findmarkers_genes, deseq2_genes)
    shared_genes <- intersect(deseq2_genes, findmarkers_genes)
    
    # Save counts to the results data frame
    results <- rbind(
      results,
      data.frame(
        Cluster = cluster_id,
        Uniquely_DESeq2 = length(unique_deseq2),
        Uniquely_FindMarkers = length(unique_findmarkers),
        Shared = length(shared_genes),
        stringsAsFactors = FALSE
      )
    )
    
    # Optionally, create and save a Venn diagram for this cluster
    venn_output <- venn.diagram(
      x = list(DESeq2 = deseq2_genes, FindMarkers = findmarkers_genes),
      filename = file.path(cluster_dir, paste0("VennDiagram_", cluster_id, ".png")),
      col = "black",
      fill = c("blue", "red"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.2,
      main = paste("Venn Diagram for", cluster_id)
    )
  } else {
    cat("Skipping cluster", cluster_id, "- missing files.\n")
  }
}

# Save the results to a CSV file
write.csv(results, file = output_file, row.names = FALSE)

cat("Results saved to", output_file, "\n")
