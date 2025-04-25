# Figure 1 plots

# This script generates:
# Panel 1H, sample tag read averages*

# * denotes that the final plot was made in Prism with output generated in R

# Load packages -----------------------------------------------------------

library(Seurat)
library(plotly)
library(tidyverse)

# Load data ---------------------------------------------------------------

load("all_10312024.RData.RData")


# Figure 1H, sample tag reads ---------------------------------------------

# Make df of sample tag reads and relevant metadata

source("Scripts/Functions/make_stdf.R")

df <- make_stdf(all)

# Splitting the dataframe by both 'st' and 'cart'
IDslist <- split(df, list(df$st, df$cart))

# Calculating the mean for each split group
mean <- sapply(IDslist, function(x) {
  # Find the names of the numeric columns
  numeric_cols <- sapply(x, is.numeric)
  
  # Ignore character columns
  numeric_cols <- numeric_cols & !sapply(x, is.character)
  
  # Find the means of the numeric columns
  means <- apply(x[, numeric_cols], 2, mean, na.rm = TRUE)
  
  # Return a named vector of means
  names(means) <- names(x)[numeric_cols]
  return(means)
})

# save the table
write.csv(mean, 'stReadsMean_bycart.csv')

# analyze ST specificity for HC samples

target_cols <- str_detect(colnames(mean), "SampleTag0[2468]_mm")
mean_hc <- mean[, target_cols]
barcode_counts <- mean_hc[c("st2","st4","st6","st8"),]

# Define the correct barcode for each sample
correct_barcodes <- c("st2", "st4", "st6", "st8")

# Add a pseudocount to avoid division by zero
pseudocount <- 1

# Initialize columns
barcode_counts$correct <- NA
barcode_counts$incorrect_mean <- NA
barcode_counts$fold_enrichment <- NA

# Loop to compute values
for (i in seq_len(nrow(barcode_counts))) {
  sample_row <- barcode_counts[i, ]
  correct_col <- correct_barcodes[i]
  
  # Extract correct and incorrect values
  correct_val <- sample_row[[correct_col]]
  incorrect_vals <- sample_row[ , setdiff(names(sample_row)[-1], c(correct_col, "Sample"))]
  incorrect_mean <- mean(as.numeric(incorrect_vals), na.rm = TRUE)
  
  # Apply pseudocount and compute fold enrichment
  barcode_counts$correct[i] <- correct_val
  barcode_counts$incorrect_mean[i] <- incorrect_mean + pseudocount
  barcode_counts$fold_enrichment[i] <- correct_val / (incorrect_mean + pseudocount)
}

# Run Wilcoxon signed-rank test against mu = 1
wilcox_result <- wilcox.test(barcode_counts$fold_enrichment, mu = 1, alternative = "greater")

# View summary and p-value
print(barcode_counts[, c("Sample", "correct", "incorrect_mean", "fold_enrichment")])
cat("Wilcoxon signed-rank test p-value:", wilcox_result$p.value, "\n")
