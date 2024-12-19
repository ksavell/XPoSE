#' Prep RRHO
#'
#' @param df_list list of pseudobulk and findmarker results combined into a df per cluster
#'
#' @return
#' @export
#'
#' @examples
prep_rrho_list <- function(df_list) {
  # Initialize a list to store the processed data frames
  processed_list <- list()
  
  # Loop through each element in the list
  for (name in names(df_list)) {
    df <- df_list[[name]]
    # Ensure the data frame is complete cases
    df <- df[complete.cases(df), ]
    
    # Calculate DDE for pseudobulk
    df[[paste0("DDE_", "pseudobulk")]] <- ifelse(
      is.na(df[["log2FoldChange"]]) | is.na(df[["padj"]]), 
      NA,  # Return NA if either log2FoldChange or padj is NA
      ifelse(df[["log2FoldChange"]] > 0, 
             -log10(df[["padj"]]), 
             ifelse(df[["log2FoldChange"]] == 0, 
                    0,
                    log10(df[["padj"]])))
    )
    
    df[[paste0("DDE_", "findmarkers")]] <- ifelse(
      is.na(df[["avg_log2FC"]]) | is.na(df[["p_val_adj"]]), 
      NA,  # Return NA if either avg_log2FC or padj is NA
      ifelse(df[["avg_log2FC"]] > 0, 
             -log10(df[["p_val_adj"]]), 
             ifelse(df[["avg_log2FC"]] == 0, 
                    0,
                    log10(df[["p_val_adj"]])))
    )
    
    # Store the processed data frame back in the list
    processed_list[[name]] <- df
  }
  
  return(processed_list)
}
