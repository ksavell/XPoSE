#' Calculate zscores from counts
#'
#' @param combined_df dataframe with the values and metadata
#' @param reference_group what is the group to z score to
#' @param val_col column of the values to z score
#' @param file_n what file to name
#'
#' @return df of combined proportions
#' @export 
#'
#' @examples
calc_zscore <- function(combined_df, reference_group = "Homecage",
                        val_col = "percent", file_n = NULL) {
        
        clust_list <- split(combined_df, combined_df[["cluster_name"]])
 
        for (i in 1:length(clust_list)) {
                current_clust <- clust_list[[i]]
                
                # Subset the data for the reference group
                reference_subset <- subset(current_clust, group == reference_group)
                
                reference_mean <- mean(reference_subset[[val_col]])
                
                reference_sd <- sd(reference_subset[[val_col]])
                
                current_clust[["z_score"]] <- (current_clust[[val_col]] - reference_mean) / reference_sd
                
                clust_list[[i]] <- current_clust
        }

# Combine the modified groups back into a single data frame
combined_df <- do.call(rbind, clust_list)

write.csv(combined_df, file = file_n)
return(combined_df)
}