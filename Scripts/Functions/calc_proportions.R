
#' Calculate proportion tables for populations
#'
#' @param seur_obj seurat object
#' @param fact1 first factor, paired with clstr_split if T
#' @param fact2 second factor
#' @param clstr_split logical to split results by cluster
#' @param file_n file name, include .csv
#'
#' @return csv file of proportions
#' @export
#'
#' @examples
#' 
calc_proportions <- function(seur_obj, fact1 = NULL, fact2 = NULL, 
                             clstr_split = FALSE, file_n = NULL){
        
        if (clstr_split == TRUE) {
                seur_obj@meta.data[["temp"]] <- paste0(seur_obj@meta.data[[fact1]],
                                                       "_", seur_obj@active.ident)
                data <- prop.table(table(seur_obj@meta.data[["temp"]], 
                                         seur_obj@meta.data[[fact2]]),
                                   margin = fact2) # this margin is broken?
        }
        else {
                data <- prop.table(table(seur_obj@meta.data[[fact1]], 
                                      seur_obj@meta.data[[fact2]]))
        }
        
        write.csv(data, file = file_n)
}

