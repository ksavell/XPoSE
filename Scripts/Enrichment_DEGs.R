# Load libraries ----------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(reshape2)

# Load data ---------------------------------------------------------------

load("glut.RData")
load("gaba.RData")

# Load functions ----------------------------------------------------------

source("Scripts/Functions/group_downsample.R")

source("Scripts/Functions/single_factor_DESeq.R")

source("Scripts/Functions/tally_iterations.R")

props <- as.numeric(c(0.035, 0.07, 0.105, 0.175, 0.245)) # 1x, 2x, 3x, 5x, 7x 

iterations <- 100

seur_obj <- gaba

seur_obj$experience <- ifelse(seur_obj$group == "Homecage", "HC", "NC")

clusters <- unique(seur_obj$cluster_name)

for (cl in clusters) {
  
  # make a folder of the cluster name
  dir_path <- paste0(cl,"/")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # no downsample
  deseq2_results <- single_factor_DESeq(object = seur_obj,
                                        comp_vect = c("experience", "NC", "HC"),
                                        cluster = cl,
                                        min_cell = 1)
  
  deseq2_results$score_column <- ifelse(deseq2_results$padj < 0.05, 1, 0)
  
  tally_file <- paste0(cl,"/experience_tally.csv")
  write.csv(deseq2_results, file = tally_file)
  
  props_list <- list()
  
  for (pr in props){
    
    all_indices <- list()
    all_seeds <- numeric(iterations)
    
    for (it in 1:iterations) {
      seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
      all_seeds[it] <- seed
      chosen_cells <- group_downsample(seur_obj, 
                                       group_to_subset = "Active", 
                                       group_to_blend = "Non-active",
                                       frac = pr, # change fraction here
                                       seed = seed)
      all_indices[[it]] <- chosen_cells
    }
    
    # Save the indices and seeds to a file
    save(all_indices, all_seeds, file = paste0(cl,"/downsampled_indices_and_seeds_",pr,".RData"))
    
    # run DESeq2 for each iteration
    results <- list()
    
    for (j in 1:iterations) {
      chosen_cells <- all_indices[[j]]
      seurat_subset <- subset(seur_obj, # change Seurat object here
                              cells = chosen_cells) 
      deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                            comp_vect = c("experience", "NC", "HC"),
                                            cluster = cl,min_cell = 1)
      results[[j]] <- deseq2_results
    }
    
    save(results, file = paste0(cl,"/results_",pr,".RData"))
    #tally_iterations(results, cl, prop_frac = pr)
    
  }
}

