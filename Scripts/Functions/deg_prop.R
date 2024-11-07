deg_prop <- function(seur_obj, iterations = 100, group_to_subset = "Active",
                     subset_frac, comp_vect, cluster){

# now 0.1 enrichment
for (i in 1:iterations) {
  seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
  all_seeds[i] <- seed
  chosen_cells <- group_downsample(seur_obj, 
                                   group_to_subset = "Active", 
                                   frac = subset_frac, # change fraction here
                                   seed = seed)
  all_indices[[i]] <- chosen_cells
}

# Save the indices and seeds to a file
save(all_indices, all_seeds, file = paste0("downsampled_indices_and_seeds_",subset_frac,".RData"))

# run DESeq2 for each iteration
results <- list()

for (i in 1:iterations) {
  chosen_cells <- all_indices[[i]]
  seurat_subset <- subset(seur_obj, # change Seurat object here
                          cells = chosen_cells) 
  deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                        comp_vect = c("experience", "NC", "HC"),
                                        cluster = "ITL23",
  )
  results[[i]] <- deseq2_results
}
file_name <- paste0(cluster,"_results_",subset_frac,".RData")
save(results, file = file_name)
tally_iterations(results, paste0(cluster,"_results_",subset_frac,"_",cluster))
}
