# Load necessary libraries
library(Seurat)
library(dplyr)
library(readr)


# Load the Seurat object from an RData file
load("/Users/holmesar/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/Robjects/all_10312024.RData")
source("Scripts/Functions/single_factor_DESeq.R")
source("Scripts/Functions/group_downsample.R")

# Define key variables
iterations <- 100
seur_obj <- all
seur_obj$experience <- ifelse(seur_obj$group == "Homecage", "HC", "NC")
clusters <- unique(seur_obj$cluster_name)
target_percentages <- c(0, 5, 10, 25, 50, 75, 100)  # Adjust percentages as needed


# Loop 1: define percent subsetting for 100 iterations --------------------


# First Loop: Calculate nuclei to keep and output RData file


# Loop through different target percentages
for (target_percentage in target_percentages) {
  
  # Initialize storage for this percentage
  all_indices <- list()
  all_seeds <- numeric(iterations)
  
  for (it in 1:iterations) {
    seed <- sample(1:10000, 1)  # Use a new seed for each iteration
    all_seeds[it] <- seed
    
    # Call group_downsample_edited
    chosen_cells <- group_downsample(
      seur_obj = seur_obj,
      group_to_subset = "Active",
      group_to_blend = "Non-active",
      percent = target_percentage,
      seed = seed
    )
    
    # Store chosen cells in the list
    all_indices[[it]] <- chosen_cells
  }
  
  # Save indices and seeds for this percentage in the working directory
  save(all_indices, all_seeds, file = paste0("downsampled_indices_and_seeds_", target_percentage, ".RData"))

}



# Loop 2: Run DESeq2 for each percent/iteration combo ---------------------


# Second Loop: Load downsampled indices and seeds, process results for each percentage, cluster, and iteration


# Loop through target percentages
for (target_percentage in target_percentages) {
  
  # Load the downsampled indices and seeds for the target percentage
  load(paste0(target_percentage, "/downsampled_indices_and_seeds_", target_percentage, ".RData"))
  
  # Loop through clusters
  for (cl in clusters) {
    
    # Create a directory for the cluster and percentage
    cluster_percentage_dir <- paste0(cl, "/", target_percentage, "/")
    if (!dir.exists(cluster_percentage_dir)) {
      dir.create(cluster_percentage_dir, recursive = TRUE)
    }
    
    # Initialize a list to store results for each iteration
    results <- list()
    
    # Loop through iterations and run DESeq2 for each
    for (it in 1:iterations) {
      
      # Extract chosen cells for the current iteration
      chosen_cells <- all_indices[[it]]
      
      # Subset the Seurat object
      seurat_subset <- subset(seur_obj, cells = chosen_cells)
      
      # Run DESeq2
      deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                            comp_vect = c("experience", "NC", "HC"),
                                            cluster = cl,
                                            min_cell = 1)
      # pull out only result tibble from the deseq2_result list
      results[[it]] <- deseq2_results$results
    }
    
    # Save the DESeq2 results for the current cluster and percentage
    save(results, file = paste0(cluster_percentage_dir, "PercentActiveInNC_experience_NC_HC_results_", target_percentage, "_percent.RData"))
  }
}
