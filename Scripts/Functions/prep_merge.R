#' Prep Merge
#' 
#' Ensures all tibbles in result list align so they can be merged. Doesn't 
#' requires sorted data so it is a bit slower.
#'
#' @param res_lst list of result tibbles
#' @param object the seurot object used to make those tibbles
#' @param gene_vect names of genes to use as rows
#'
#' @return the corrected result list
#' @export
#'
#' @examples
prep_merge <- function(res_lst, object, gene_vect = object@assays[["RNA"]]@counts@Dimnames[[1]]){
    timer <- proc.time()
    
    #will be returned at function's end
    return_lst <- list()
    
    ind <- 0
    
    for (cluster in names(res_lst)) {
        #Tells user where function is in correction process
        ind <- ind + 1
        print(paste("On Index ", ind, "/", length(res_lst), sep = ""), quote = FALSE)
        
        #allows us to take rows from the base tibble
        cls_frm <- data.frame(res_lst[[cluster]])
        
        #makes what will be new tibble and fills genes
        new_tbl <- data.frame(matrix(nrow = length(gene_vect),
                                     ncol = ncol(res_lst[[cluster]])))
        
        #Adds some necessary data to new tibble
        new_tbl[, 1] <- gene_vect
        colnames(new_tbl) <- colnames(cls_frm)
        rownames(cls_frm) <- cls_frm[, "gene"]
        
        for (i in which(new_tbl[, 1] %in% res_lst[[cluster]]$gene)) {
            new_tbl[i, ] <- cls_frm[new_tbl[i, "gene"], ]
        }
        new_tbl[, 1] <- gene_vect
        
        #returns to tibble and apppends to list
        new_tbl <- tibble(new_tbl)
        
        return_lst <- list.append(return_lst, new_tbl)
    }
    
    #ensures names are same
    names(return_lst) <- names(res_lst)
    
    #tells user how long function took
    elpd <- (proc.time() - timer)[[3]]
    print(paste("Process took", round(elpd, digits = 3), "seconds"), 
          quote = FALSE)
    
    #returns new list
    return(return_lst)
}
