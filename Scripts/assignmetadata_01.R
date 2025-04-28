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

# Load data ---------------------------------------------------------------

# Read gene count tables
counts_c1 <- read.csv("~/input/Combined_cart1-good-only_RSEC_MolsPerCell.csv", 
                      skip = 7, row.names = 1)
counts_c2 <- read.csv("~/input/Combined_cart2-good-only_RSEC_MolsPerCell.csv", 
                      skip = 7, row.names = 1)

# Read sample tag calls
tags_c1 <- read.csv("~/input/cart1-good-only_Sample_Tag_Calls.csv", 
                    skip = 7, row.names = 1)
tags_c2 <- read.csv("~/input/cart2-good-only_Sample_Tag_Calls.csv", 
                    skip = 7, row.names = 1)

#Read Sample tag reads
streads_c1 <- read.csv("~/input/cart1-good-only_Sample_Tag_ReadsPerCell.csv", 
                       skip = 7, row.names = 1)
streads_c2 <- read.csv("~/input/cart2-good-only_Sample_Tag_ReadsPerCell.csv", 
                       skip = 7, row.names = 1)

# Reorder all input files -------------------------------------------------

# current objects in environment to order
objects <- ls()

#order them
invisible(lapply(objects, function(name) {
  obj <- get(name, envir = .GlobalEnv)
  if (is.data.frame(obj)) {
    sorted_obj <- obj[order(rownames(obj)), ]
    assign(name, sorted_obj, envir = .GlobalEnv)
  }
}))

# Create intermediate object and bind sample tag reads and calls ----------

# create lists
streads <- list(c1 = streads_c1, c2 = streads_c2)
counts <- list(c1 = counts_c1, c2 = counts_c2)
tags <- list(t1 = subset(tags_c1, select=-c(2)),
             t2 = subset(tags_c2, select=-c(2)))

source("~/XPoSE/Scripts/Functions/create_seur.R")

# creates separate seurat object for each count table present in environment
# In this experiment, we used 2 cartridges and used sample tags # 2-9.
for (i in 1:length(streads)) {
  assign(paste("data","_c",i,sep=""),
         create_seur(counts[i], i, tags[[i]],
                     addSTreads = TRUE, streads = streads[[i]], STrange = c(2:9)),
         envir = .GlobalEnv)
  
}

# Assign experimental metadata based on cartridge and ST assignment -------

# Creating a list of both Seurat Objects for easier reference
seur_list <- list(data_c1 = data_c1, data_c2 = data_c2)

# input list of ratIDs used
ratIDs <- c(NA, "HC-3", "NC-4", "HC-2", "NC-2", "HC-4",
            "NC-3", "HC-1", "NC-1")

# input sex of each rat
rat_sex <- c(NA, "Male", "Male", "Female", "Female", "Male",
             "Male", "Female", "Female")

# creates sample tag list
st_list <- c(NA)
for(i in 2:length(ratIDs)){
  st_list <- c(st_list,paste("SampleTag0",i,"_mm",sep=""))
}

source("~/XPoSE/Scripts/Functions/assign_combine.R") # read function for assignment

for (l in 1:length(seur_list)) {
  seur_list <- assign_combine(seur_list, l, ratIDs, st_list, rat_sex)
}

combined <- merge(seur_list[[1]],seur_list[[2]])

# include ST calls that in experiment
combined <- subset(combined, (Sample_tag=="SampleTag02_mm" | 
                                 Sample_tag=="SampleTag03_mm" |
                                 Sample_tag=="SampleTag04_mm" |
                                 Sample_tag=="SampleTag05_mm" |
                                Sample_tag=="SampleTag06_mm" |
                                Sample_tag=="SampleTag07_mm" |
                                Sample_tag=="SampleTag08_mm" |
                                Sample_tag=="SampleTag09_mm"))

save(combined, file = "combined08032023.RData")


# Prep for MapMyCells -----------------------------------------------------

library(Seurat)
library(reticulate)
py_require(c("anndata"))

load("combined08032023.RData")

# Python modules
anndata <- import("anndata")
scipy <- import("scipy.sparse")
pd <- import("pandas")

obj <- combined

# pull counts and format
counts <- t(as.matrix(obj[["RNA"]]@counts))  # Make genes = columns
sparse_counts_py <- scipy$csr_matrix(counts)

genes   <- colnames(counts)
samples <- rownames(counts)

var_df <- pd$DataFrame(dict(genes = genes), index = genes)
obs_df <- pd$DataFrame(dict(samples = samples), index = samples)

# convert and save as .h5ad
countAD <- anndata$AnnData(X = sparse_counts_py, var = var_df, obs = obs_df)

countAD$write("combinedcounts.h5ad")