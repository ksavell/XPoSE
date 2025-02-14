# Load libraries ----------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(reshape2)

# Load data ---------------------------------------------------------------

load("all_10312024.RData")

# Load functions ----------------------------------------------------------

source("Scripts/Functions/group_downsample.R")

source("Scripts/Functions/single_factor_DESeq.R")

source("Scripts/Functions/tally_iterations.R")


# Define key variables ----------------------------------------------------

iterations <- 100
seur_obj <- all 
seur_obj$experience <- ifelse(seur_obj$group == "Homecage", "HC", "NC")
clusters <- unique(seur_obj$cluster_name)
ratios <- c(0.5, 1, 2)  # Define the ratios for downsampling


# loop through clusters ---------------------------------------------------


for (cl in clusters) {
  
  # Create a folder for the cluster
  dir_path <- paste0(cl, "/")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Run DESeq2 without downsampling (baseline)
  deseq2_results <- single_factor_DESeq(object = seur_obj,
                                        comp_vect = c("experience", "NC", "HC"),
                                        cluster = cl,
                                        min_cell = 1)
  deseq2_results$score_column <- ifelse(deseq2_results$padj < 0.05, 1, 0)
  
  # Save baseline results
  write.csv(deseq2_results, file = paste0(dir_path, "experience_tally.csv"))
  
  # Loop through different ratios
  for (ratio in ratios) {
    
    all_indices <- list()
    all_seeds <- numeric(iterations)
    
    for (it in 1:iterations) {
      seed <- sample(1:10000, 1)  # Use a new seed for each iteration
      all_seeds[it] <- seed
      
      # Get downsampled cells
      chosen_cells <- group_downsample_ratio(seur_obj,
                                             group_to_subset = "Stimulated",
                                             group_to_blend = "Control",
                                             ratio = ratio,
                                             seed = seed)
      all_indices[[it]] <- chosen_cells
    }
    
    # Save indices and seeds for this ratio
    save(all_indices, all_seeds, file = paste0(dir_path, "downsampled_indices_and_seeds_", ratio, ".RData"))
    
    # Run DESeq2 for each downsampled iteration
    results <- list()
    
    for (j in 1:iterations) {
      chosen_cells <- all_indices[[j]]
      seurat_subset <- subset(seur_obj, cells = chosen_cells) 
      
      deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                            comp_vect = c("experience", "NC", "HC"),
                                            cluster = cl,
                                            min_cell = 1)
      results[[j]] <- deseq2_results
    }
    
    # Save DESeq2 results for this ratio
    save(results, file = paste0(dir_path, "results_", ratio, ".RData"))
    
    # Optional: Call a tally function if needed
    # tally_iterations(results, cl, prop_frac = ratio)
  }
}