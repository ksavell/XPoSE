# degs between populations

load("~/rstudio_projects/XPoSE/all_10312024.RData")
source("Scripts/Functions/single_factor_DESeq.R")

# now run DEGs based on populations

# loop through clusters for the 3 comparisons

clusters <- unique(all$cluster_name)

for (cl in clusters) {
  tryCatch({
    # Run DESeq2 for Non-active vs Homecage
    deseq2_results <- single_factor_DESeq(object = all,
                                          comp_vect = c("group", "Non-active", "Homecage"),
                                          cluster = cl,
                                          min_cell = 1)
    deseq2_results$score_column <- ifelse(deseq2_results$padj < 0.05, 1, 0)
    deg_file <- paste0(cl, "_group_Nonactive_Homecage.csv")
    write.csv(deseq2_results, file = deg_file)
    
    # Run DESeq2 for Active vs Homecage
    deseq2_results <- single_factor_DESeq(object = all,
                                          comp_vect = c("group", "Active", "Homecage"),
                                          cluster = cl,
                                          min_cell = 1)
    deg_file <- paste0(cl, "_group_Active_Homecage.csv")
    write.csv(deseq2_results, file = deg_file)
    
    # Run DESeq2 for Active vs Non-active
    deseq2_results <- single_factor_DESeq(object = all,
                                          comp_vect = c("group", "Active", "Non-active"),
                                          cluster = cl,
                                          min_cell = 1)
    deg_file <- paste0(cl, "_group_Active_Nonactive.csv")
    write.csv(deseq2_results, file = deg_file)
  }, error = function(e) {
    # Print a message and continue to the next cluster
    message(paste("Error for cluster:", cl, "-", e$message))
  })
}

# Process CSV files for a given comparison

# List of comparisons
comparisons <- c("group_Nonactive_Homecage", "group_Active_Homecage", "group_Active_Nonactive")

# Function to process CSV files for a given comparison
tally_differential_expression <- function(comparison) {
  # List all files matching the comparison
  csv_files <- list.files(pattern = paste0("_", comparison, ".csv$"))
  
  # Initialize tally data frame
  tally_df <- data.frame(
    Cluster = character(),
    Upregulated = integer(),
    Downregulated = integer(),
    Total = integer(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each file
  for (file in csv_files) {
    # Read the CSV
    data <- read.csv(file)
    
    # Tally upregulated and downregulated genes
    upregulated <- sum(data$log2FoldChange > 0 & data$padj < 0.05, na.rm = TRUE)
    downregulated <- sum(data$log2FoldChange < 0 & data$padj < 0.05, na.rm = TRUE)
    
    # Extract cluster name from file name
    cluster <- gsub(paste0("_", comparison, ".csv$"), "", file)
    
    # Add to tally data frame
    tally_df <- rbind(tally_df, data.frame(
      Cluster = cluster,
      Upregulated = upregulated,
      Downregulated = downregulated,
      Total = upregulated + downregulated,
      stringsAsFactors = FALSE
    ))
  }
  
  # Save the tally as a CSV
  output_file <- paste0("Tally_", comparison, ".csv")
  write.csv(tally_df, file = output_file, row.names = FALSE)
  
  # Return the tally data frame
  return(tally_df)
}

# Loop over comparisons
for (comparison in comparisons) {
  message(paste("Processing comparison:", comparison))
  tally_results <- tally_differential_expression(comparison)
  print(tally_results)  # Optional: Print results for each comparison
}
