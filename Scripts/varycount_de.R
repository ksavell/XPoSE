# F5 count DEG

load("~/Projects/XPoSE/all_10312024.RData")
source("Scripts/Functions/group_downsample_counts.R")
source("Scripts/Functions/single_factor_DESeq.R")

library(tidyverse)

# vary the absolute number of clusters and run DEG
# 532 nuclei at baseline

# Define key variables
iterations <- 100
seur_obj <- all
clusters <- unique(seur_obj$cluster_name)
target_numbers <- c(5, 10, 25, 50, 100, 200, 300, 400, 500)  # Adjust numbers as needed


# Loop 1: define number subsetting for 100 iterations --------------------


# First Loop: Calculate nuclei to keep and output RData file

for (cl in clusters) {
# Loop through different target numbers
  cluster_number_dir <- file.path(cl)
  if (!dir.exists(cluster_number_dir)) {
    dir.create(cluster_number_dir, recursive = TRUE)
  }
for (target_number in target_numbers) {
  
  # Initialize storage for this number
  all_indices <- list()
  all_seeds <- numeric(iterations)
  
  for (it in 1:iterations) {
    seed <- sample(1:10000, 1)  # Use a new seed for each iteration
    all_seeds[it] <- seed
    
    # Call group_downsample_counts
    chosen_cells <- group_downsample_counts(
      seur_obj = seur_obj,
      group_to_subset = "Active",
      group_to_blend = "Non-active",
      number = target_number,
      seed = seed,
      cluster = cl
    )
    
    # Store chosen cells in the list
    all_indices[[it]] <- chosen_cells
  }
  
  # Save indices and seeds for this number in the working directory
  save(all_indices, all_seeds, file = paste0(cl,"/", cl, "_downsampled_indices_and_seeds_count_", target_number, ".RData"))
    
    # Initialize a list to store results for each iteration
    results <- list()
    
    # Loop through iterations and run DESeq2 for each
    for (j in seq_along(all_indices)) {
      
      # Extract chosen cells for the current iteration
      chosen_cells <- all_indices[[j]]
      
      result <- tryCatch({
        seurat_subset <- subset(seur_obj, cells = chosen_cells)
        deseq2_results <- single_factor_DESeq(
          object = seurat_subset,
          comp_vect = c("group", "Active", "Non-active"),
          cluster = cl,
          min_cell = 0
        )
        deseq2_results$deseq_results  # Extract result tibble
      }, error = function(e) {
        message(paste("Error in cluster", cl, "at iteration", j, "for", target_number, "count:", e$message))
        NULL  # Return NULL for failed iterations
      })
      
      results[[j]] <- result
    }
    
    # Save the DESeq2 results for the current cluster and number
    save(results, file = file.path(cl, paste0("/Number_group_Active_Nonactive_results_", target_number, "_count.RData")))
  }
}

