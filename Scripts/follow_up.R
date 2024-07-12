# Additional analysis

# consulted with MLG to talk through suggestions

# # Run DEGs with 100%, 50%, or 25% of the Active population

# suggestion through reading:
  # calculate delta-variance for each biological replicate through the Libra package

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

# prep metadata  -------------------------------------------------------

# make experience metadata

glut$experience <- ifelse(glut$group == "Homecage", "HC",
                          "NC")


# Experience, no downsample -----------------------------------------------

deseq2_results <- single_factor_DESeq(object = glut,
                                      comp_vect = c("experience", "NC", "HC"),
                                      cluster = "ITL23"
)

deseq2_results$score_column <- ifelse(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0, 1,
                                      ifelse(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < 0, -1, 0))

# 0.035 prop ----------------------------------------------------------------

iterations <- 100
all_indices <- list()
all_seeds <- numeric(iterations)

for (i in 1:iterations) {
  seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
  all_seeds[i] <- seed
  chosen_cells <- group_downsample(glut, 
                                   group_to_subset = "Active", 
                                   frac = 0.035, # change fraction here
                                   seed = seed)
  all_indices[[i]] <- chosen_cells
}

# Save the indices and seeds to a file
save(all_indices, all_seeds, file = "downsampled_indices_and_seeds_0p035.RData")

# run DESeq2 for each iteration
itl23_results <- list()

for (i in 1:iterations) {
  chosen_cells <- all_indices[[i]]
  seurat_subset <- subset(glut, # change Seurat object here
                          cells = chosen_cells) 
  deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                        comp_vect = c("experience", "NC", "HC"),
                                        cluster = "ITL23",
                                        )
  itl23_results[[i]] <- deseq2_results
}
save(itl23_results, file = "itl23_results_0p035.RData")
tally_iterations(itl23_results, "ITL23", prop_frac = "p_035")


# 0.1 ---------------------------------------------------------------------

# now 0.1 enrichment
for (i in 1:iterations) {
  seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
  all_seeds[i] <- seed
  chosen_cells <- group_downsample(glut, 
                                   group_to_subset = "Active", 
                                   frac = 0.1, # change fraction here
                                   seed = seed)
  all_indices[[i]] <- chosen_cells
}

# Save the indices and seeds to a file
save(all_indices, all_seeds, file = "downsampled_indices_and_seeds_0p1.RData")

# run DESeq2 for each iteration
itl23_results <- list()

for (i in 1:iterations) {
  chosen_cells <- all_indices[[i]]
  seurat_subset <- subset(glut, # change Seurat object here
                          cells = chosen_cells) 
  deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                        comp_vect = c("experience", "NC", "HC"),
                                        cluster = "ITL23",
  )
  itl23_results[[i]] <- deseq2_results
}
save(itl23_results, file = "itl23_results_0p1.RData")
tally_iterations(itl23_results, "ITL23", prop_frac = "p_10")

#deg_prop(glut,subset_frac = 0.1, comp_vect = c("experience","NC","HC"), cluster = "ITL23")


# 0.25 prop ----------------------------------------------------------------

# now 0.25 enrichment
for (i in 1:iterations) {
  seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
  all_seeds[i] <- seed
  chosen_cells <- group_downsample(glut, 
                                   group_to_subset = "Active", 
                                   frac = 0.25, # change fraction here
                                   seed = seed)
  all_indices[[i]] <- chosen_cells
}

# Save the indices and seeds to a file
save(all_indices, all_seeds, file = "downsampled_indices_and_seeds_0p25.RData")

# run DESeq2 for each iteration
itl23_results <- list()

for (i in 1:iterations) {
  chosen_cells <- all_indices[[i]]
  seurat_subset <- subset(glut, # change Seurat object here
                          cells = chosen_cells) 
  deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                        comp_vect = c("experience", "NC", "HC"),
                                        cluster = "ITL23",
  )
  itl23_results[[i]] <- deseq2_results
}
save(itl23_results, file = "itl23_results_0p25.RData")
tally_iterations(itl23_results, "ITL23", prop_frac = "p_25")


# 0.5 prop ----------------------------------------------------------------

for (i in 1:iterations) {
  seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
  all_seeds[i] <- seed
  chosen_cells <- group_downsample(glut, 
                                   group_to_subset = "Active", 
                                   frac = 0.5, # change fraction here
                                   seed = seed)
  all_indices[[i]] <- chosen_cells
}

# Save the indices and seeds to a file
save(all_indices, all_seeds, file = "downsampled_indices_and_seeds_0p5.RData")

# run DESeq2 for each iteration
itl23_results <- list()

for (i in 1:iterations) {
  chosen_cells <- all_indices[[i]]
  seurat_subset <- subset(glut, # change Seurat object here
                          cells = chosen_cells) 
  deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                        comp_vect = c("experience", "NC", "HC"),
                                        cluster = "ITL23",
  )
  itl23_results[[i]] <- deseq2_results
}
save(itl23_results, file = "itl23_results_0p5.RData")
tally_iterations(itl23_results, "ITL23", prop_frac = "p_50")

# 0.75 prop ----------------------------------------------------------------

for (i in 1:iterations) {
  seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
  all_seeds[i] <- seed
  chosen_cells <- group_downsample(glut, 
                                   group_to_subset = "Active", 
                                   frac = 0.75, # change fraction here
                                   seed = seed)
  all_indices[[i]] <- chosen_cells
}

# Save the indices and seeds to a file
save(all_indices, all_seeds, file = "downsampled_indices_and_seeds_0p75.RData")

# run DESeq2 for each iteration
itl23_results <- list()

for (i in 1:iterations) {
  chosen_cells <- all_indices[[i]]
  seurat_subset <- subset(glut, # change Seurat object here
                          cells = chosen_cells) 
  deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                        comp_vect = c("experience", "NC", "HC"),
                                        cluster = "ITL23",
  )
  itl23_results[[i]] <- deseq2_results
}
save(itl23_results, file = "itl23_results_0p75.RData")
tally_iterations(itl23_results, "ITL23", prop_frac = "p_75")