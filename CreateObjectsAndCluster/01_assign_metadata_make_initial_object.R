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
counts_c1 <- read.table("~/MakeObject/Combined_cart1-good-only_RSEC_MolsPerCell.csv", 
                        sep = ",", skip = 7, header = TRUE, row.names = 1)
counts_c2 <- read.table("~/MakeObject/Combined_cart2-good-only_RSEC_MolsPerCell.csv", 
                        sep = ",", skip = 7, header = TRUE, row.names = 1)

# Read sample tag calls
tags_c1 <- read.table("~/MakeObject/cart1-good-only_Sample_Tag_Calls.csv", 
                      sep = ",", skip = 7, header = TRUE, row.names = 1)
tags_c2 <- read.table("~/MakeObject/cart2-good-only_Sample_Tag_Calls.csv", 
                      sep = ",", skip = 7, header = TRUE, row.names = 1)

#Read Sample tag reads
streads_c1 <- read.table("~/MakeObject/cart1-good-only_Sample_Tag_ReadsPerCell.csv", 
                         sep = ",", skip = 7, header = TRUE, row.names = 1)
streads_c2 <- read.table("~/MakeObject/cart2-good-only_Sample_Tag_ReadsPerCell.csv", 
                         sep = ",", skip = 7, header = TRUE, row.names = 1)


# Reorder all input files -------------------------------------------------

# function to order rows (nuclei barcodes)
order_rows <- function(data) {
  ordered_data <- data[order(row.names(data)), ]
  return(ordered_data)
}

counts_c1 <- order_rows(counts_c1)
counts_c2 <- order_rows(counts_c2)
tags_c1 <- order_rows(tags_c1)
tags_c2 <- order_rows(tags_c2)
streads_c1 <- order_rows(streads_c1)
streads_c2 <- order_rows(streads_c2)

# Create intermediate object and bind sample tag reads and calls ----------
# In this experiment, we used 2 cartridges and used sample tags #2-9.

# create Seurat objects for each cartridge
# transpose the counts since SB output is opposite of Seurat input
data_c1 <- CreateSeuratObject(counts = t(counts_c1), project = "C1")
data_c2 <- CreateSeuratObject(counts = t(counts_c2), project = "C2")

# add sample tag calls as metadata
data_c1$Sample_tag <- cbind(tags_c1$Sample_Tag)
data_c2$Sample_tag <- cbind(tags_c2$Sample_Tag)

# Define the Seurat objects
seurat_objects <- list(data_c1, data_c2)

# Define the SampleTag numbers used in the experiment
sample_tags <- paste0("ST", 2:9)

# Loop over the Seurat objects and SampleTags

