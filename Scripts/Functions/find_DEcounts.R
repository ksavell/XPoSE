find_DEcounts <- function(directory, cluster, de_path = "_group_Active_Homecage", feature_list, control_suffix = "HC") {
  
  # Construct file paths
  
  dds_file <- file.path(directory, paste0(cluster, de_path, ".RDS"))
  coexp_file <- file.path(directory, paste0(cluster, de_path, ".csv"))
  
  # Load DESeq2 dataset
  if (!file.exists(dds_file) | !file.exists(coexp_file)) {
    stop("Missing input files for cluster: ", cluster)
  }
  dds <- readRDS(dds_file)
  coexp <- read.csv(coexp_file, row.names = 1)
  
  # Extract normalized counts
  count <- counts(dds, normalized = TRUE)
  
  # Define relevant variables
  FCvar <- "log2FoldChange"
  PVvar <- "padj"
  
  # Filter by feature list
  filtered_count <- as.data.frame(count[rownames(count) %in% feature_list, ])
  
  # Identify control samples
  control_samples <- colnames(count)[grepl(control_suffix, colnames(count))]
  
  # Calculate log2 fold change relative to control group
  control_avg <- rowMeans(filtered_count[, control_samples, drop = FALSE], na.rm = TRUE)
  
  for (sample in colnames(filtered_count)) {
    new_column_name <- paste0("log2FC_", sample)
    filtered_count[[new_column_name]] <- log2(filtered_count[[sample]] / control_avg)
  }
  
  # Copy adjusted p-values from coexp table
  idx <- match(rownames(filtered_count), rownames(coexp))
  filtered_count$adjpval <- coexp[[PVvar]][idx]
  filtered_count$adjpvalT <- -log10(filtered_count$adjpval)
  
  # Save as CSV
  output_file <- file.path(directory, paste0("Log2FC_Values_", cluster, de_path, ".csv"))
  write.csv(filtered_count, file = output_file, row.names = TRUE)
  
  message("Saved file: ", output_file)
}
