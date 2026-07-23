# Figure 2 scripts

# generates data for figure 2 panels:
# Panel F: nuclei counts per sample
# Panel G: Mean reads by XPoSE-tag and associated stats
# Cell type verification by MapMyCells annotation

# Info --------------------------------------------------------------------

# This script generates data for panels F-H

# Loading -----------------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)

source("Scripts/Functions/make_stdf.R")
source("Scripts/Functions/calc_prop.R")

# load in unclustered object that is output of createobject_01.R

# F2F ---------------------------------------------------------------------

rat_orig_table <- table(obj$ratID, obj$orig.ident)
rat_orig_df <- as.data.frame(rat_orig_table)

# save
write.csv(rat_orig_df, "f2f_ratID_vs_origident_table.csv", row.names = FALSE)

# F2G ---------------------------------------------------------------------

## means
df <- make_stdf(obj)

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

# save
write.csv(mean, 'f2g_stReadsMean_bycart.csv')

## stats 

# list to hold results
wilcox_results <- list()

# Define  samples and corresponding incorrect tags
samples <- list(
  "SampleTag08" = c("SampleTag02_reads", "SampleTag04_reads", "SampleTag06_reads"),
  "SampleTag04" = c("SampleTag02_reads", "SampleTag06_reads", "SampleTag08_reads"),
  "SampleTag02" = c("SampleTag04_reads", "SampleTag06_reads", "SampleTag08_reads"),
  "SampleTag06" = c("SampleTag02_reads", "SampleTag04_reads", "SampleTag08_reads")
)

# loop
for (sample_base in names(samples)) {
  
  assigned_tag <- paste0(sample_base, "_mm")   # Assigned Sample_tag
  correct_reads_var <- paste0(sample_base, "_reads")  # Correct reads column
  
  # Subset cells
  cells <- WhichCells(obj, expression = Sample_tag == assigned_tag)
  
  # Fetch correct and incorrect tag reads
  correct <- FetchData(obj, vars = correct_reads_var)[cells, 1]
  incorrect <- rowMeans(FetchData(obj, vars = samples[[sample_base]])[cells, ])
  
  # Run Wilcoxon test
  test <- wilcox.test(correct, incorrect, paired = TRUE, alternative = "greater")
  
  # Save the results
  wilcox_results[[sample_base]] <- data.frame(
    sample = sample_base,
    p_value = test$p.value,
    statistic = test$statistic,
    method = test$method,
    alternative = test$alternative
  )
}

# Combine all results into a single data frame
wilcox_summary <- do.call(rbind, wilcox_results)

# Save
write.csv(wilcox_summary, "f2g_stats.csv", row.names = FALSE)

# F2H ---------------------------------------------------------------------

obj_celltype <- calc_prop(seur_obj = obj, 
                          fact1 = 'ratID',
                          fact2 = 'celltype',
                          fact3 = 'orig.ident')

write.csv(obj_celltype, file = "f2h_celltype_by_rat_cart.csv")
