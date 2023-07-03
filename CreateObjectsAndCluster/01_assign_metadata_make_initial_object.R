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
library(tidyverse)

# Load data ---------------------------------------------------------------

# set working directory where the raw data files are found
setwd("~/MakeObject")

# Read gene count tables
counts.C1 <- read.table("~/MakeObject/Combined_cart1-good-only_RSEC_MolsPerCell.csv", 
                        sep = ",", skip = 7, header = TRUE, row.names = 1)
counts.C2 <- read.table("~/MakeObject/Combined_cart2-good-only_RSEC_MolsPerCell.csv", 
                        sep = ",", skip = 7, header = TRUE, row.names = 1)

# Read sample tag calls
tags.C1 <- read.table("~/MakeObject/cart1-good-only_Sample_Tag_Calls.csv", 
                      sep = ",", skip = 7, header = TRUE, row.names = 1)
tags.C2 <- read.table("~/MakeObject/cart2-good-only_Sample_Tag_Calls.csv", 
                      sep = ",", skip = 7, header = TRUE, row.names = 1)

#Read Sample tag reads
STreads.C1 <- read.table("~/MakeObject/cart1-good-only_Sample_Tag_ReadsPerCell.csv", 
                         sep = ",", skip = 7, header = TRUE, row.names = 1)
STreads.C2 <- read.table("~/MakeObject/cart2-good-only_Sample_Tag_ReadsPerCell.csv", 
                         sep = ",", skip = 7, header = TRUE, row.names = 1)

#reorder so binding works
tags.C1 <- tags.C1[order(row.names(tags.C1)),]
tags.C2 <- tags.C2[order(row.names(tags.C2)),]

counts.C1 <- counts.C1[order(row.names(counts.C1)),]
counts.C2 <- counts.C2[order(row.names(counts.C2)),]


# Create intermediate object and bind sample tag reads and calls ----------
# This is custom since we used 2 cartridges and used sample tags #2-9.

# create Seurat objects for each cartridge
# transpose the counts since SB output is opposite of Seurat input
data.C1 <- CreateSeuratObject(counts = t(counts.C1), project = "C1")
data.C2 <- CreateSeuratObject(counts = t(counts.C2), project = "C2")

# add sample tag calls as metadata
data.C1$Sample_tag <- cbind(tags.C1$Sample_Tag)
data.C2$Sample_tag <- cbind(tags.C2$Sample_Tag)

#Add Sample_tag reads as metadata to Seurat object
data.C1$ST2_reads <- cbind(STreads.C1$SampleTag02_mm.stAbO)
data.C1$ST3_reads <- cbind(STreads.C1$SampleTag03_mm.stAbO)
data.C1$ST4_reads <- cbind(STreads.C1$SampleTag04_mm.stAbO)
data.C1$ST5_reads <- cbind(STreads.C1$SampleTag05_mm.stAbO)
data.C1$ST6_reads <- cbind(STreads.C1$SampleTag06_mm.stAbO)
data.C1$ST7_reads <- cbind(STreads.C1$SampleTag07_mm.stAbO)
data.C1$ST8_reads <- cbind(STreads.C1$SampleTag08_mm.stAbO)
data.C1$ST9_reads <- cbind(STreads.C1$SampleTag09_mm.stAbO)

data.C2$ST2_reads <- cbind(STreads.C2$SampleTag02_mm.stAbO)
data.C2$ST3_reads <- cbind(STreads.C2$SampleTag03_mm.stAbO)
data.C2$ST4_reads <- cbind(STreads.C2$SampleTag04_mm.stAbO)
data.C2$ST5_reads <- cbind(STreads.C2$SampleTag05_mm.stAbO)
data.C2$ST6_reads <- cbind(STreads.C2$SampleTag06_mm.stAbO)
data.C2$ST7_reads <- cbind(STreads.C2$SampleTag07_mm.stAbO)
data.C2$ST8_reads <- cbind(STreads.C2$SampleTag08_mm.stAbO)
data.C2$ST9_reads <- cbind(STreads.C2$SampleTag09_mm.stAbO)

