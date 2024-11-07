# F5

source("Scripts/Functions/single_factor_DESeq.R")

clusters <- unique(all$cluster_name)

pairs_list <- list(
  c("group", "Non-active", "Homecage"),
  c("group", "Active", "Homecage"),
  c("group", "Active", "Non-active")
)

# Initialize the results data frame
results <- data.frame(Category = character(), Observation = character(), Value = numeric(), stringsAsFactors = FALSE)

for (pair in pairs_list) {
  for (cl in clusters) {
    # Use tryCatch to handle cases where single_factor_DESeq() fails or returns NULL
    deseq2_results <- tryCatch({
      single_factor_DESeq(object = all, comp_vect = pair, cluster = cl, min_cell = 1)
    }, error = function(e) {
      NULL  # Return NULL if an error occurs
    })
    
    # Skip this iteration if deseq2_results is NULL or empty
    if (is.null(deseq2_results) || nrow(deseq2_results) == 0) {
      next  # Skip to the next iteration
    }
    
    # Check if necessary columns exist in deseq2_results
    if (!("padj" %in% colnames(deseq2_results)) || !("log2FoldChange" %in% colnames(deseq2_results))) {
      warning(paste("Columns 'padj' or 'log2FoldChange' missing for cluster:", cl))
      next  # Skip to the next iteration
    }
    
    # Ensure padj and log2FoldChange are numeric
    deseq2_results$padj <- as.numeric(deseq2_results$padj)
    deseq2_results$log2FoldChange <- as.numeric(deseq2_results$log2FoldChange)
    
    # Add up_score and dn_score columns to the DESeq2 results
    deseq2_results$up_score <- ifelse(!is.na(deseq2_results$padj) & deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0, 1, 0)
    deseq2_results$dn_score <- ifelse(!is.na(deseq2_results$padj) & deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < 0, -1, 0)
    
    # Append the up_score and dn_score to the results data frame
    results <- rbind(results, data.frame(Category = paste(pair, collapse = "_"), Observation = cl, Value = sum(deseq2_results$up_score, na.rm = TRUE)))
    results <- rbind(results, data.frame(Category = paste(pair, collapse = "_"), Observation = cl, Value = sum(deseq2_results$dn_score, na.rm = TRUE)))
  }
}

# Check the results
print(results)

row_colors <- c(
  'CTL6' = '#2C8CB9',
  'PTL5' = '#0A5B8C',
  'ITL23' = '#41B75F',
  'ITL5' = '#5DBFC1', 
  'ITL6' = '#3A8F87',
  'NPL56' = '#3C9E64',
  'CTL6b' = '#6F499D',
  'Pvalb' = '#E66027',
  'Sst' = '#F8991D',
  'Meis2' = '#C52126',
  'Vip' = '#A669AB',
  'Lamp5' = '#DB808C',
  'SstChodl' = '#B0B235',
  'PvalbChand' = '#AD6C49')

ggplot(results, aes(x = Category, y = Value, color = Observation)) +
  geom_point(size = 3, position = position_jitter(width = 0.2, height = 0)) +  # Adds slight horizontal jitter for better visibility
  scale_color_manual(values = row_colors) +  # Generates unique colors for each 'Observation'
  labs(
       x = "Category",
       y = "Score") +
  theme_minimal() +
  theme(
    legend.position = "none",  # Removes the legend for a cleaner plot
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotates x-axis labels for readability
  )


