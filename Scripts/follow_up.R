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

source("Scripts/Functions/group_downsample.R")

# downsample group:active to 50% 

glut_downsampled <- group_downsample(seur_obj = glut, 
                                     group_to_subset = "Active", 
                                     frac = 0.5, 
                                     bio_rep = "ratID")

source("Scripts/Functions/single_factor_DESeq.R")

# Run DESeq2 on the downsampled object

# need a loop for this for cluster
ITL23 <- single_factor_DESeq(object = glut_downsampled,
                                        comp_vect = c("group", "Active", "Non-active"),
                                        cluster = "ITL23")

# suggestion for iteration

# Run 100 iterations and save the indices of the cells chosen
iterations <- 100
all_indices <- list()

for (i in 1:iterations) {
  seed <- i # You can change the seed if needed
  chosen_cells <- group_downsample(glut, 
                                   group_to_subset = "Active", 
                                   frac = 0.5, 
                                   bio_rep = "ratID", 
                                   seed = 222)
  all_indices[[i]] <- chosen_cells
}

# Save the indices to a file if needed
save(all_indices, file = "downsampled_indices.RData")

# If you need to use any specific iteration for further analysis
iteration_to_use <- 1
chosen_cells <- all_indices[[iteration_to_use]]
seurat_subset <- subset(seurat_object, cells = chosen_cells)