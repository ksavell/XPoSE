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
calc_prop <- function(seur_obj, fact1 = 'ratID', fact2 = 'cluster_name',
                      fact3 = NULL, file_n = NULL) {
        
        if(is.null(fact3)){
                count_df <- data.frame(table(seur_obj@meta.data[[fact1]],
                                         seur_obj@meta.data[[fact2]]))
                colnames(count_df) <- c(fact1, fact2, "count")
        }
        else {
                count_df <- data.frame(table(seur_obj@meta.data[[fact1]],
                                         seur_obj@meta.data[[fact2]],
                                         seur_obj@meta.data[[fact3]]))
                colnames(count_df) <- c(fact1, fact2, fact3, "count")
        }
        
        count_df[[paste0(fact1,"_",fact3)]] <- paste0(count_df[[fact1]], count_df[[fact3]])
        
IDslist <- split(count_df, count_df[[paste0(fact1,"_",fact3)]])

for (i in 1:length(IDslist)) {
        current_group <- IDslist[[i]]
        
        total_count <- sum(current_group[["count"]])
        current_group[["percent"]] <- current_group[["count"]] / total_count
        
        IDslist[[i]] <- current_group
}

# Combine the modified groups back into a single data frame
combined_df <- do.call(rbind, IDslist)
combined_df <- combined_df[complete.cases(combined_df[c("count", "percent")]), ]

write.csv(combined_df, file = file_n)
return(combined_df)
}