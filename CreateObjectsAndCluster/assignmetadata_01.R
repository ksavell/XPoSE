# Assign metadata and create combined initial seurat object
# Authors:  Katherine Savell and Padmashri Saravanan

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

combined <- merge(data_c1, data_c2) # combine data from both cartridges

# use assign_combine.R for any other metadata assignments