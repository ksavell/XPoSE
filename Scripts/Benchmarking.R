# Benchmarking

# Load libraries ----------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)


# Load data ---------------------------------------------------------------

load("glut.RData")
load("gaba.RData")

# Animal representation in the dataset

# One advantage of multiplexing is you can observe 
# animal overrepresentation within a pooled capture

# Percent representation of each animal in the dataset

source("Scripts/Functions/calc_prop.R")

glut_group <- calc_prop(seur_obj = glut, fact1 = 'ratID',
                        fact2 = "orig.ident",
                        file_n = "glut_group.csv")


# correct vs not for animal representation --------------------------------

table(glut$ratID,glut$cluster_name)

# Assuming 'glut' is your Seurat object

# Create an empty list to store the subset Seurat objects
subset_seurat_objs <- list()

# Set the number of iterations
iterations <- 250

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


