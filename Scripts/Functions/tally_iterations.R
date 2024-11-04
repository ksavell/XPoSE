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
tally_iterations <- function(results_list, j) {

# Initialize a list to store DESeq2 results with scores for each iteration
scored_results <- list()

# Process each iteration's DESeq2 results
for (i in 1:100) {
  # Get the DESeq2 results for this iteration
  deseq2_df <- as.data.frame(results_list[[j]][["results"]][[i]])
  
  # Add the gene names as a column for merging purposes
  rownames(deseq2_df) <- deseq2_df$gene
  
  # Add a score column based on the conditions
  score_column <- ifelse(deseq2_df$padj < 0.05, 1, 0)
  
  # Add the score column to the data frame with a unique name for the iteration
  col_name <- paste0("score_iteration_", i)
  deseq2_df[[col_name]] <- score_column
  
  # Store the scored results in the list, keeping only gene and the score column
  scored_results[[i]] <- deseq2_df[, c("gene", col_name)]
}
print(paste("Processed", length(scored_results), "iterations"))

# Merge the results into a single data frame using the gene names as the key
merged_results_df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), scored_results)

# Restore the gene names as rownames and remove the temporary 'gene' column
rownames(merged_results_df) <- merged_results_df$gene 
merged_results_df <- merged_results_df[, -which(names(merged_results_df) == "gene")]

return(merged_results_df)
}


