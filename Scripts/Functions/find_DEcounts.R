find_DEcounts <- function(dds, coexp, cluster, feature_list, factor_metadata, factor_levels){
  
  source("Scripts/Functions/sort_df.R")
  
  # Extract all normalized counts for each sample
  count <- counts(dds, normalized=TRUE)
  
  # Define the cluster and genes to pull (sig. upregulated)
  FCvar <- "log2FoldChange"
  PVvar <- "padj"
  
  # Filter by feature list
  filtered_count <- as.data.frame(count[rownames(count) %in% feature_list, ])
  
  # Calculate averages based on the provided factor levels
  averages <- sapply(factor_levels, function(level) {
    sample_cols <- factor_metadata$sample[factor_metadata$factor == level]
    rowMeans(filtered_count[, sample_cols, drop = FALSE], na.rm = TRUE)
  })
  
  # Add the averages as new columns
  filtered_count <- cbind(filtered_count, averages)
  
  # Calculate fold changes for each sample relative to the factor-based average
  for (i in seq_along(colnames(filtered_count[, 1:ncol(count)]))) {
    new_column_name <- paste0("FC_", colnames(filtered_count)[i])
    filtered_count[new_column_name] <- log2(filtered_count[, i] / filtered_count[, "average"])
  }
  
  # Copy over the adjusted p-value from the coexp table
  idx <- match(rownames(filtered_count), rownames(coexp))
  filtered_count$adjpval <- coexp[[PVvar]][idx]
  filtered_count$adjpvalT <- -log10(filtered_count$adjpval)
  
  # Save the filtered count as a CSV
  write.csv(filtered_count, file = paste0("ExampleFeatures_", cluster, ".csv"))
}
