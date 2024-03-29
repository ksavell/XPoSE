#' Prep Merge (Fast)
#' 
#' Ensures all tibbles in result list align so they can be merged. Requires 
#' sorted data as it uses a binary search algorithm.
#'
#' @param res_lst a list of result tibbles
#' @param object Seurat object analysis was performed on
#' @param gene_vect names of genes to use as rows
#'
#' @return the corrected result list
#' @export
#'
#' @examples
#' #Prepares tibbles in res object.
#' res <- prep_merge(res, gs)
prep_merge_fast <- function(res_lst, object, gene_vect = object@assays[["RNA"]]@counts@Dimnames[[1]]){

    source("~/Documents/GitHub/XPoSE/Scripts/Functions/bi_search.R")
    #Tracks time function took to run
    timer <- proc.time()
    #will be returned at function's end
    return_lst <- list()
    
    ind <- 0
    
    #goes through entire results list
    for (cluster in names(res_lst)) {
        #Tells user where function is in correction process
        ind <- ind + 1
        print(paste("On Index ", ind, "/", length(res_lst), sep = ""), quote = FALSE)
        
        #allows us to take rows from the base tibble
        cls_frm <- data.frame(res_lst[[cluster]])
        
        #new tibble that will added to list
        new_tbl <- data.frame(matrix(nrow = length(gene_vect),
                                     ncol = ncol(res_lst[[cluster]])))
        
        #Adds necessary data to new tibble
        new_tbl[, 1] <- gene_vect
        colnames(new_tbl) <- colnames(cls_frm)
        rownames(cls_frm) <- cls_frm[, "gene"]
        
        #fills shared rows w/ base data
        for (i in 1:nrow(new_tbl)) {
            if(!is.na(bi_search(cls_frm[, 1], new_tbl[i, 1]))){
                #replaces row
                new_tbl[i, ] <- cls_frm[bi_search(cls_frm[, 1], new_tbl[i, 1]), ]
            }
        }
        
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
