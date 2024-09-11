#' Calculate proportions from counts
#'
#' @param seur_obj seurat object
#' @param fact1 the first factor to consider, usually ratID
#' @param fact2 the second factor to consider, like group, sex, or cartridge
#'
#' @return df of combined proportions
#' @export 
#'
#' @examples
calc_prop <- function(seur_obj, fact1 = 'ratID', fact2 = 'celltype', fact3 = NULL) {
  
  if (is.null(fact3)) {
    # Case when fact3 is not provided, only use fact1 and fact2
    count_df <- data.frame(table(seur_obj@meta.data[[fact1]], seur_obj@meta.data[[fact2]]))
    colnames(count_df) <- c(fact1, fact2, "count")
  } else {
    # Case when fact3 is provided, use fact1, fact2, and fact3
    count_df <- data.frame(table(seur_obj@meta.data[[fact1]], seur_obj@meta.data[[fact2]], seur_obj@meta.data[[fact3]]))
    colnames(count_df) <- c(fact1, fact2, fact3, "count")
  }
  
  # Now split the data by the combination of fact1 and fact3 (if fact3 exists)
  if (!is.null(fact3)) {
    count_df[[paste0(fact1, "_", fact3)]] <- paste0(count_df[[fact1]], "_", count_df[[fact3]])
    IDslist <- split(count_df, count_df[[paste0(fact1, "_", fact3)]])
  } else {
    IDslist <- split(count_df, count_df[[fact1]])
  }
  
  # Iterate through each group, calculate the percentage
  for (i in 1:length(IDslist)) {
    current_group <- IDslist[[i]]
    
    total_count <- sum(current_group[["count"]])
    current_group[["percent"]] <- current_group[["count"]] / total_count
    
    IDslist[[i]] <- current_group
  }
  
  # Combine the modified groups back into a single data frame
  combined_df <- do.call(rbind, IDslist)
  combined_df <- combined_df[complete.cases(combined_df[c("count", "percent")]), ]
  
  return(combined_df)
}

