#' Make Metadata
#' 
#' Swipes the rat:factor pairs from the big table in order to make the metadata
#'
#' @param tbl_lst list of the psuedobulk counts
#' @param object Seurat object to perform analysis on
#'
#' @return returns metadata as data frame
#' @export
#'
#' @examples
#' #Make metadata for the glut object
#' meta_data <- prep_DESeq(my_tbls, glut)
prep_DESeq <- function(tbl_lst, object){
    
    #splits names into list of the ratID's and groups
    temp <- sapply(names(tbl_lst[[1]]), strsplit, ":")
    justF <- c()
    
    #Fills justF w/ only the comparisons in the split
    for (i in 1:length(temp)){
        justF <- append(justF, temp[[i]][2]) 
    }
    
    #finds factor used in the list for user
    for (name in colnames(object@meta.data)){
        if (justF[1] %in% unique(object@meta.data[name][, 1])){
            factor <- name
        }
    }
    
    #Makes data frame and labels it 
    met_fram <- data.frame(justF)
    rownames(met_fram) <- names(tbl_lst[[1]])
    colnames(met_fram) <- factor  #<- Will likely need to be changed if multiple
                                  #   factors wanted
    
    #checks that the metadata is set up correctly
    cat("There should be at least one 'TRUE' under this line:\n")
    cat(as.character(unique(met_fram[,1]) %in% unique(object@meta.data[factor][, 1])), "\n")
    cat("\n")
    
    #returns frame
    return(met_fram)
}
