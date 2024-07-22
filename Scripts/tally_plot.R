# read in tally csv's

library(ggplot2)
library(reshape2)
library(dplyr)

source("Scripts/Functions/tally_iterations.R")

# load results first 

clstr <- "ITL23"

files <- list.files(clstr, pattern = "^results_.*\\.RData$", full.names = TRUE)
data_list <- lapply(files, function(file) {
  env <- new.env()
  load(file, envir = env)
  as.list(env)
})

names(data_list) <- basename(files)

# remove the file extension from the names
names(data_list) <- sub("\\.RData$", "", names(data_list))

# now tally the results
for (j in names((data_list))) {
  merged_results_df <- tally_iterations(data_list, j)
  write.csv(merged_results_df, file = paste0(clstr,"/",j,"_iteration_tally.csv"))
}











# count occurrence of +1 in each row for each tally and write into a new df

count_observed <- data.frame(p035 = (rowSums(ITL23p_035iteration_tally == 1, na.rm = TRUE)),
                             p10 = (rowSums(ITL23p_10iteration_tally == 1, na.rm = TRUE)),
                             p25 = (rowSums(ITL23p_25iteration_tally == 1, na.rm = TRUE)),
                             p50 = (rowSums(ITL23p_50iteration_tally == 1, na.rm = TRUE)),
                             p75 = (rowSums(ITL23p_75iteration_tally == 1, na.rm = TRUE)),
                             p100 = (ITL23_100iteration_tally$score_column),
                             row.names = rownames(ITL23p_035iteration_tally))

# subset rows that have at least one +1 in the 100 iteration tally
row_keep <- which(count_observed$p100 == 1)

count_observed <- count_observed[row_keep,]

# Function to calculate mean and SEM
calc_mean_sem_count <- function(column) {
  
  # Calculate mean
  mean_value <- mean(column)
  
  # Calculate SEM
  sem_value <- sd(column) / sqrt(length(column))
  
  return(c(mean = mean_value, SEM = sem_value))
}

# Apply the function to each column in the data frame
result <- lapply(count_observed, calc_mean_sem_count)

# Convert the result to a data frame
result_df <- as.data.frame(do.call(rbind, result))

write_csv(result_df, "ITL23_prop_observed.csv")
