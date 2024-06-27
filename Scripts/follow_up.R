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


# Animal representation in dataset ----------------------------------------

# One advantage of multiplexing is you can observe 
# animal overrepresentation within a pooled capture

# Percent representation of each animal in the dataset

source("Scripts/Functions/calc_prop.R")

glut_group <- calc_prop(seur_obj = glut, fact1 = 'ratID',
                        fact2 = "orig.ident",
                        file_n = "glut_group.csv")


# correct vs not for animal representation --------------------------------

table(glut$ratID,glut$cluster_name)

# Create an empty list to store the subset Seurat objects
subset_seurat_objs <- list()

# Set the number of iterations
iterations <- 100

# Set the subset size
subset_size <- 100

# Perform iterations
for (i in 1:iterations) {
  # Sample indices with replacement
  subset_indices <- sample(1:ncol(glut), size = subset_size, replace = TRUE)
  # Subset the Seurat object
  subset_seurat_obj <- glut[, subset_indices]
  # Append the subset object to the list
  subset_seurat_objs[[i]] <- subset_seurat_obj
}



# downsample active -------------------------------------------------------

table(glut$group, glut$ratID)

# run the non-downsampled results

deseq2_results <- single_factor_DESeq(object = glut,
                                      comp_vect = c("group", "Active", "Non-active"),
                                      cluster = "ITL23",
)

deseq2_results$score_column <- ifelse(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0, 1,
                                      ifelse(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < 0, -1, 0))

# now time to downsample

source("Scripts/Functions/group_downsample.R")

source("Scripts/Functions/single_factor_DESeq.R")

source("Scripts/Functions/tally_iterations.R")

iterations <- 100
all_indices <- list()
all_seeds <- numeric(iterations)

# ITL23 as example, 0.5 ---------------------------------------------------

for (i in 1:iterations) {
  seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
  all_seeds[i] <- seed
  chosen_cells <- group_downsample(glut, 
                                   group_to_subset = "Active", 
                                   frac = 0.5, # change fraction here
                                   bio_rep = "ratID", 
                                   seed = seed)
  all_indices[[i]] <- chosen_cells
}

# Save the indices and seeds to a file
save(all_indices, all_seeds, file = "downsampled_indices_and_seeds_ITL23_0-5.RData")

itl23_results <- list()

for (i in 1:iterations) {
  chosen_cells <- all_indices[[i]]
  seurat_subset <- subset(glut, # change Seurat object here
                          cells = chosen_cells) 
  deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                        comp_vect = c("group", "Active", "Non-active"),
                                        cluster = "ITL23",
                                        )
  itl23_results[[i]] <- deseq2_results
}

tally_iterations(itl23_results, "ITL23_0-5_")

# ITL23 as example, 0.25 ---------------------------------------------------

for (i in 1:iterations) {
  seed <- sample(1:10000, 1)  # Use a random seed for each iteration to ensure distinct subsets
  all_seeds[i] <- seed
  chosen_cells <- group_downsample(glut, 
                                   group_to_subset = "Active", 
                                   frac = 0.25, # change fraction here
                                   bio_rep = "ratID", 
                                   seed = seed)
  all_indices[[i]] <- chosen_cells
}

# Save the indices and seeds to a file
save(all_indices, all_seeds, file = "downsampled_indices_and_seeds_ITL23_0-25.RData")

itl23_results <- list()

for (i in 1:iterations) {
  chosen_cells <- all_indices[[i]]
  seurat_subset <- subset(glut, # change Seurat object here
                          cells = chosen_cells) 
  deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                        comp_vect = c("group", "Active", "Non-active"),
                                        cluster = "ITL23",
  )
  itl23_results[[i]] <- deseq2_results
}

tally_iterations(itl23_results, "ITL23_0-25_")
