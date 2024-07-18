# read in tally csv's

library(ggplot2)
library(reshape2)
library(dplyr)

source("Scripts/Functions/tally_iterations.R")

# load results first 
files <- list.files("ITL23", pattern = "^results_.*\\.RData$", full.names = TRUE)
data_list <- lapply(files, function(file) {
  env <- new.env()
  load(file, envir = env)
  as.list(env)
})

# now tally the results
for (j in 1:length(data_list)) {
  results <- data_list[[j]]
  tally_iterations(results)
  write.csv(merged_results_df, file = paste0(i,"/iteration_tally.csv"))
}

tally_iterations(results, i, prop_frac = "p_035")

ITL23p_75iteration_tally <- read.csv("ITL23p_75iteration_tally.csv", row.names = 1, header = TRUE)
ITL23p_50iteration_tally <- read.csv("ITL23p_50iteration_tally.csv", row.names = 1, header = TRUE)
ITL23p_25iteration_tally <- read.csv("ITL23p_25iteration_tally.csv", row.names = 1, header = TRUE)
ITL23p_10iteration_tally <- read.csv("ITL23p_10iteration_tally.csv", row.names = 1, header = TRUE)
ITL23p_035iteration_tally <- read.csv("ITL23p_035iteration_tally.csv", row.names = 1, header = TRUE)
ITL23_100iteration_tally <- read.csv("ITL23_experience_tally.csv", row.names = 1, header = TRUE)

# order ITL23_100iteration_tally by gene column name
ITL23_100iteration_tally <- ITL23_100iteration_tally[order(ITL23_100iteration_tally$gene),]


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
