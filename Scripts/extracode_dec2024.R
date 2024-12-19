# drop rat extras

# tally the results (collapses into one value) ----------------------------


# Initialize an empty list to store results
tally_results <- list()

# Loop through directories for each excluded rat
for (excluded_rat in all_rats) {
  
  # Path to the excluded rat's directory
  rat_dir <- paste0("Excluded_", excluded_rat)
  
  # Get a list of CSV files in the directory
  csv_files <- list.files(rat_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize a data frame to store tallies for the current rat
  rat_tally <- data.frame(Cluster = character(), ExcludedRat = character(), SignificantGenes = integer())
  
  # Loop through CSV files and process each
  for (file in csv_files) {
    # Read the CSV file
    deseq_results <- read.csv(file)
    
    # Tally significant genes
    significant_genes <- sum(deseq_results$padj < 0.05, na.rm = TRUE)
    
    # Extract cluster name from file name
    cluster_name <- gsub(".*Cluster_(.+)_Excluded_.*\\.csv$", "\\1", basename(file))
    
    # Add results to the tally data frame
    rat_tally <- rbind(rat_tally, data.frame(
      Cluster = cluster_name,
      ExcludedRat = excluded_rat,
      SignificantGenes = significant_genes
    ))
  }
  
  # Append the current rat's tally to the results list
  tally_results[[excluded_rat]] <- rat_tally
}

# Combine all results into a single data frame
final_tally <- do.call(rbind, tally_results)

# Save the final tally as a CSV
write.csv(final_tally, "Final_Tally_Significant_Genes.csv", row.names = FALSE)

# Display the final tally
print(final_tally)


# plots the tally ---------------------------------------------------------

# plot it
library(ggplot2)

# Define custom colors for the excluded rats
color_mapping <- c(
  "HC-1" = "gray20", "HC-2" = "gray40", "HC-3" = "gray60", "HC-4" = "gray80",
  "NC-1" = "orange", "NC-2" = "darkorange", "NC-3" = "orangered", "NC-4" = "darkorange3"
)

# Create the plot
ggplot(final_tally, aes(x = Cluster, y = SignificantGenes, fill = ExcludedRat)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = color_mapping) +
  labs(
    title = "Number of Differentially Expressed Genes by Cluster",
    x = "Cluster",
    y = "Number of Significant Genes",
    fill = "Excluded Rat"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


