# read in tally csv's

library(ggplot2)
library(reshape2)


ITL23p_75iteration_tally <- read.csv("ITL23p_75iteration_tally.csv", row.names = 1, header = TRUE)
ITL23p_50iteration_tally <- read.csv("ITL23p_50iteration_tally.csv", row.names = 1, header = TRUE)
ITL23p_25iteration_tally <- read.csv("ITL23p_25iteration_tally.csv", row.names = 1, header = TRUE)
ITL23p_10iteration_tally <- read.csv("ITL23p_10iteration_tally.csv", row.names = 1, header = TRUE)
ITL23p_035iteration_tally <- read.csv("ITL23p_035iteration_tally.csv", row.names = 1, header = TRUE)

# count occurance of +1 and -1 in each row for each tally and write into a new df

count_observed <- data.frame(p035 = (rowSums(ITL23p_035iteration_tally == 1, na.rm = TRUE) + 
                               rowSums(ITL23p_035iteration_tally == -1, na.rm = TRUE)),
                             p10 = (rowSums(ITL23p_10iteration_tally == 1, na.rm = TRUE) + 
                               rowSums(ITL23p_10iteration_tally == -1, na.rm = TRUE)),
                             p25 = (rowSums(ITL23p_25iteration_tally == 1, na.rm = TRUE) + 
                               rowSums(ITL23p_25iteration_tally == -1, na.rm = TRUE)),
                             p50 = (rowSums(ITL23p_50iteration_tally == 1, na.rm = TRUE) + 
                               rowSums(ITL23p_50iteration_tally == -1, na.rm = TRUE)),
                             p75 = (rowSums(ITL23p_75iteration_tally == 1, na.rm = TRUE) + 
                               rowSums(ITL23p_75iteration_tally == -1, na.rm = TRUE)))

count_observed <- as.data.frame(lapply(count_observed, function(x) x / 100))

prop_observed <- data.frame(p035 = count_observed$p035/deseq2_results$score_column,
                            p10 = count_observed$p10/deseq2_results$score_column,
                            p25 = count_observed$p25/deseq2_results$score_column,
                            p50 = count_observed$p50/deseq2_results$score_column,
                            p75 = count_observed$p75/deseq2_results$score_column)

#broken
for (i in 1:ncol(count_observed)) {
  prop_observed[, i] <- ifelse(count_observed[, i] == 0, 0, count_observed[, i] / deseq2_results$score_column)
}

