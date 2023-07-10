# Assign metadata and create combined initial seurat object
# Authors: Katherine Savell and Padmashri Saravanan

# Info --------------------------------------------------------------------

# This script takes 3 input files per cartridge from the alignment pipeline
# in SB
# count table RSEC_MolsPerCell.csv
# sample tag reads Sample_Tag_ReadsPerCell.csv
# sample tag calls Sample_Tag_Calls.csv

# Load packages -----------------------------------------------------------
library(Seurat)
#library(tidyverse)

# Load data ---------------------------------------------------------------

# set working directory where the raw data files are found
setwd("~/Desktop/Analysis/XPoSE_seq_cartridge")

# Read gene count tables
counts_c1 <- read.csv("~/Desktop/Analysis/XPoSE_seq_cartridge/Combined_cart1-good-only_RSEC_MolsPerCell.csv", 
                        skip = 7, row.names = 1)
counts_c2 <- read.csv("~/Desktop/Analysis/XPoSE_seq_cartridge/Combined_cart2-good-only_RSEC_MolsPerCell.csv", 
                        skip = 7, row.names = 1)

# Read sample tag calls
tags_c1 <- read.csv("~/Desktop/Analysis/XPoSE_seq_cartridge/cart1-good-only_Sample_Tag_Calls.csv", 
                    skip = 7, row.names = 1)
tags_c2 <- read.csv("~/Desktop/Analysis/XPoSE_seq_cartridge/cart2-good-only_Sample_Tag_Calls.csv", 
                    skip = 7, row.names = 1)

#Read Sample tag reads
streads_c1 <- read.csv("~/Desktop/Analysis/XPoSE_seq_cartridge/cart1-good-only_Sample_Tag_ReadsPerCell.csv", 
                       skip = 7, row.names = 1)
streads_c2 <- read.csv("~/Desktop/Analysis/XPoSE_seq_cartridge/cart2-good-only_Sample_Tag_ReadsPerCell.csv", 
                       skip = 7, row.names = 1)


# Reorder all input files -------------------------------------------------

# current objects in environment to order
objects <- ls()

# function to order by rows or specified column
sort_dataframe <- function(df, col_sort = NULL) {
        if (is.null(col_sort)) {
                # Sort based on row names
                sorted_df <- df[order(row.names(df)), ]
        } else {
                if (col_sort %in% colnames(df)) {
                        # Sort based on the specified column
                        sorted_df <- df[order(df[, col_sort]), ]
                } else {
                        stop("Column name not found!")
                }
        }
        
        return(sorted_df)
}

#order them
invisible(lapply(objects, function(name) {
        df <- get(name)  # Get the data frame from the environment
        modified_df <- sort_dataframe(df)  # Apply the function to the data frame
        assign(name, modified_df, envir = .GlobalEnv)  # Update the data frame in the environment
}))

# Create intermediate object and bind sample tag reads and calls ----------
# In this experiment, we used 2 cartridges and used sample tags #2-9.

# create Seurat objects for each cartridge
# transpose the counts since SB output is opposite of Seurat input
data_c1 <- CreateSeuratObject(counts = t(counts_c1), project = "c1")
data_c2 <- CreateSeuratObject(counts = t(counts_c2), project = "c2")

# add sample tag calls as metadata
data_c1$stag <- cbind(tags_c1$Sample_Tag)
data_c2$stag <- cbind(tags_c2$Sample_Tag)

add_STreads <- function(seurat_obj, reads_obj, STrange = c(1:12)) {
  for (i in STrange) {
    colName <- paste0("SampleTag", i, "_reads")
    seurat_obj[[colName]] <- cbind(reads_obj[[paste0("SampleTag", sprintf("%02d", i), "_mm.stAbO")]])
  }
  return(seurat_obj)
}

streads <- list(c1 = streads_c1, c2 = streads_c2)
counts <- list(c1 = counts_c1, c2 = counts_c2)
tags <- list(t1 = subset(tags_c1, select=-c(2)),
             t2 = subset(tags_c2, select=-c(2)))

source("~/Documents/GitHub/XPoSE/Scripts/Functions/create_seur.R")
# creates separate seurat object for each count table present in environment
for (i in 1:length(streads)) {
  assign(paste("data","_c",i,sep=""),
         create_seur(counts[i], i, tags[[i]],
                     addSTreads = TRUE, streads = streads[[i]], STrange = c(2:9)),
         envir = .GlobalEnv)
  
}


# Add Sample_tag reads as metadata to Seurat objects
data_c1 <- add_STreads(data_c1, streads_c1, STrange =  2:9)
data_c2 <- add_STreads(data_c2, streads_c2, STrange =  2:9)

combined <- merge(data_c1, data_c2)

# BROKEN Assign metadata ---------------------------------------------------------
# Create a cart_stag metadata
combined$cart_stag <- paste0(combined$orig.ident, combined$stag, sep = "_")
Idents(combined) <- 'cart_stag'

# Might need to exclude multiplet and undetermined
combined_filt <- subset(combined, idents = c("c1SampleTag02_mm_", "c2SampleTag02_mm_",
                                             "c1SampleTag03_mm_", "c2SampleTag03_mm_",
                                             "c1SampleTag04_mm_", "c2SampleTag04_mm_",
                                             "c1SampleTag05_mm_", "c2SampleTag05_mm_",
                                             "c1SampleTag06_mm_", "c2SampleTag06_mm_",
                                             "c1SampleTag07_mm_", "c2SampleTag07_mm_",
                                             "c1SampleTag08_mm_", "c2SampleTag08_mm_"))

# add_metadata <- function(seurat_obj, metadata_df, split_col = cart_stag, factor1 = NULL, 
#                          factor2 = NULL, factor3 = NULL, factor4 = NULL, factor5 = NULL) {
#   # Check if the column to split by exists in the Seurat object
#   if (!(split_col %in% colnames(seurat_obj@meta.data))) {
#     stop(paste("Column", split_col, "not found in seurat_obj"))
#   }
#   
#   # Get unique values of the split column
#   unique_vals <- unique(seurat_obj@meta.data[[split_col]])
#   
#   # Make a list to store split objects
#   split_data <- vector("list", length(unique_vals))
#   
#   # Iterate over unique values and create split objects
#   for (i in seq_along(unique_vals)) {
#     # Filter rows based on the current unique value
#     split_rows <- seurat_obj@meta.data[[split_col]] == unique_vals[i]
#     
#     # Create a new Seurat object for the split rows
#     split_obj <- subset(seurat_obj, subset = split_rows)
#     
#     # Extract metadata for the split object
#     metadata <- metadata_df[rownames(metadata_df) == unique_vals[i], ]
#   
#     # Assign metadata to the split object
#     split_obj$factor1 <- metadata$factor1
#     split_obj$factor2 <- metadata$factor2
#     split_obj$factor3 <- metadata$factor3
#     split_obj$factor4 <- metadata$factor4
#     split_obj$factor5 <- metadata$factor5
#     
#     # Update the split data
#     split_data[[i]] <- split_obj
#   }
#   
# # Recombine the split objects into a single Seurat object
# seurat_obj <- do.call(Seurat, split_data)
#   
#   # Return the combined Seurat object
#   return(seurat_obj)
# }
# 
# combined_meta <- add_metadata(combined_filt, split_col = 'cart_stag', metadata = metadata, 
#                               factor1 = 'group', factor2 = 'ratID')

# Let's try again
# This is not working yet but I think it's promising
split_add_merge <- function(seurat_obj, meta_name, new_meta_df){
  # Check if the column to split by exists in the Seurat object
     if (!(meta_name %in% colnames(seurat_obj@meta.data))) {
       stop(paste("Column", meta_name, "not found in seurat_obj"))
     }
  # Get the unique metadata values
  unique_meta <- unique(seurat_obj@meta.data[[meta_name]])
  
  # Initialize an empty list to store split Seurat objects
  seurat_list <- list()
  
  # Iterate over the unique metadata values
  for (i in seq_along(unique_meta)) {
    # Split the Seurat object by metadata
    temp_obj <- subset(seurat_obj, subset = seurat_obj@meta.data[[meta_name]] == unique_meta[i])
    # Add new metadata to the split object
    for (j in 1:ncol(new_meta_df)) {
      temp_obj[[colnames(new_meta_df)[j]]] <- new_meta_df[i, j]
    }
    # Add the updated Seurat object to the list
    seurat_list[[i]] <- temp_obj
  }
  
  # Recombine the Seurat objects
  new_seurat_obj <- seurat_list[[1]]
  for (i in seq_along(seurat_list)[-1]) {
    new_seurat_obj <- merge(new_seurat_obj, y = seurat_list[[i]])
  }
  
  # Return the new Seurat object
  return(new_seurat_obj)
}

combined_meta <- split_add_merge(combined_filt, meta_name = 'cart_stag', new_meta_df = metadata)
