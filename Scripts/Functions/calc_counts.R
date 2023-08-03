#' Calculate count tables for populations
#'
#' @param seur_obj seurat object
#' @param fact1 first factor, rows
#' @param fact2 second factor, columns
#' @param fact3 third factor, final splitter
#' @param # clstr_split logical to split fact1 by cluster
#' @param file_n file name, include .csv
#'
#' @return csv file of counts
#' @export
#'
#' @examples
#' 
calc_counts <- function(seur_obj, fact1 = NULL, fact2 = NULL,
                        fact3 = NULL, file_n = NULL){
        
        #if (clstr_split == TRUE) {
                #seur_obj@meta.data[["temp"]] <- paste0(seur_obj@meta.data[[fact1]],
                #                                      "_", seur_obj@active.ident)
                #data <- table(seur_obj@meta.data[["temp"]], 
                #                        seur_obj@meta.data[[fact2]])
        #}
       # else {
        
        if(is.null(fact3)){
                data <- data.frame(table(seur_obj@meta.data[[fact1]],
                                         seur_obj@meta.data[[fact2]]))
                colnames(data) <- c(fact1, fact2, "count")
        }
        else {
                data <- data.frame(table(seur_obj@meta.data[[fact1]],
                                         seur_obj@meta.data[[fact2]],
                                         seur_obj@meta.data[[fact3]]))
                colnames(data) <- c(fact1, fact2, fact3, "count")
        }
        
        #}
        
        write.csv(data, file = file_n)
        return(data)
}

