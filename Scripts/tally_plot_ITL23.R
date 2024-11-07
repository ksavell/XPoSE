# read in tally csv's

library(ggplot2)
library(reshape2)
library(dplyr)

source("Scripts/Functions/tally_iterations.R")
load("all_10312024.RData")

#clusters <- unique(all$cluster_name)

# loop through clusters

#for (cl in clusters) {
#clstr <- cl 

files <- list.files(pattern = "ITL23_results_.*\\.RData$", full.names = TRUE)
data_list <- lapply(files, function(file) {
  env <- new.env()
  load(file, envir = env)
  as.list(env)
})

names(data_list) <- basename(files)

# remove the file extension from the names
names(data_list) <- sub("\\.RData$", "", names(data_list))

# now tally the results and save that csv 
for (j in names(data_list)) {
  merged_results_df <- tally_iterations(data_list, j)
  write.csv(merged_results_df, 
            file = paste0(j, "_iteration_tally.csv"))
}



