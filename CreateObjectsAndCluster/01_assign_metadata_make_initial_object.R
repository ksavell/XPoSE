# Assign metadata and create combined initial seurat object
# Authors: Katherine Savell

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
setwd("~/MakeObject")

# Read gene count tables
counts_c1 <- read.csv("~/MakeObject/Combined_cart1-good-only_RSEC_MolsPerCell.csv", 
                        skip = 7, row.names = 1)
counts_c2 <- read.csv("~/MakeObject/Combined_cart2-good-only_RSEC_MolsPerCell.csv", 
                        skip = 7, row.names = 1)

# Read sample tag calls
tags_c1 <- read.csv("~/MakeObject/cart1-good-only_Sample_Tag_Calls.csv", 
                    skip = 7, row.names = 1)
tags_c2 <- read.csv("~/MakeObject/cart2-good-only_Sample_Tag_Calls.csv", 
                    skip = 7, row.names = 1)

#Read Sample tag reads
streads_c1 <- read.csv("~/MakeObject/cart1-good-only_Sample_Tag_ReadsPerCell.csv", 
                       skip = 7, row.names = 1)
streads_c2 <- read.csv("~/MakeObject/cart2-good-only_Sample_Tag_ReadsPerCell.csv", 
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
data_c1$sample_tag <- cbind(tags_c1$Sample_Tag)
data_c2$sample_tag <- cbind(tags_c2$Sample_Tag)

addSTReads <- function(seuratObj, readsObj, tagPrefix = 'ST', STrange = c(1:12)) {
  for (i in STrange) {
    colName <- paste0(tagPrefix, i, "_reads")
    seuratObj[[colName]] <- cbind(readsObj[[paste0("SampleTag", sprintf("%02d", i), "_mm.stAbO")]])
  }
  return(seuratObj)
}

# Add Sample_tag reads as metadata to Seurat objects
data_c1 <- addSTReads(data_c1, streads_c1, STrange =  2:9)
data_c2 <- addSTReads(data_c2, streads_c2, STrange =  2:9)


