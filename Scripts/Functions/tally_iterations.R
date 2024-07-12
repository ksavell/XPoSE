# tally deg iterations

#' Title
#'
#' @param results_list list of deseq2 results from each iteration
#' @param cluster name of the cluster
#'
#' @return
#' @export
#'
#' @examples
tally_iterations <- function(results_list, cluster, prop_frac) {

# Initialize a list to store DESeq2 results with scores for each iteration
scored_results <- list()

# Process each iteration's DESeq2 results
for (i in 1:length(results_list)) {
  # Get the DESeq2 results for this iteration
  deseq2_df <- as.data.frame(results_list[[i]])
  
  # Add the gene names as a column for merging purposes
  rownames(deseq2_df) <- deseq2_df$gene
  
  # Add a score column based on the conditions
  score_column <- ifelse(deseq2_df$padj < 0.05 & deseq2_df$log2FoldChange > 0, 1,
                         ifelse(deseq2_df$padj < 0.05 & deseq2_df$log2FoldChange < 0, -1, 0))
  
  # Add the score column to the data frame with a unique name for the iteration
  col_name <- paste0("score_iteration_", i)
  deseq2_df[[col_name]] <- score_column
  
  # Store the scored results in the list
  scored_results[[i]] <- deseq2_df[, c("gene", col_name)]
}

# Merge the results into a single data frame using the gene names as the key
merged_results_df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), scored_results)

# Restore the gene names as rownames and remove the temporary 'gene' column
rownames(merged_results_df) <- merged_results_df$gene
merged_results_df <- merged_results_df[, -which(names(merged_results_df) == "gene")]

# save the merged results to a file

write.csv(merged_results_df, file = paste0(cluster,prop_frac,"iteration_tally.csv"))

# make a table that counts 1, 0, and -1 values for each iteration
# this is not working yet...

# # Initialize a data frame to store the counts
# summary <- data.frame(count_minus1 = 0,
#                             count_0 = 0,
#                             count_1 = 0)
# 
# # Iterate over each column and compute the frequency table
# for (col in names(merged_results_df)) {  
#   # Compute frequency table for the current column
#   col_table <- as.data.frame(table(merged_results_df[[col]]))
#   
#   # Update summary table for the current column
#   summary[col, "count_minus1"] <- col_table[col_table$Var1 == -1, "Freq"]
#   summary[col, "count_0"] <- col_table[col_table$Var1 == 0, "Freq"]
#   summary[col, "count_1"] <- col_table[col_table$Var1 == 1, "Freq"]
# }
# 
# summary <- summary[-1,]
# 
# # Print the consolidated summary table
# write.csv(summary, file = paste0(prop_frac,cluster,"iteration_summary.csv"))

}
