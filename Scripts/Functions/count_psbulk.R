#' Count Psuedobulk
#' 
#' #Makes a list containing the count tables for each of the user's chosen 
#' clusters. Can "deal" with multiple factors
#'
#' @param data_lst a list of Seurat objects with user chosen clusters
#' @param object Seurat object to perform analysis on
#'
#' @return a list of the psuedobulk counts
#' @export
#'
#' @examples
#' #Make list of count tables for glut object
#' my_tbls <- count_pseudobulk(clust_lst, glut)
count_psbulk <- function(data_lst, object, incl_all = FALSE){
    
    #Will hold all the tables
    furniture <- list(matrix(1:length(data_lst[[1]]@assays[["RNA"]]@counts@Dimnames[[1]]),
                             ncol = 1))
    
    #finds factor used in the list for user
    for (name in colnames(object@meta.data)){
        if (unique(data_lst[[1]]@active.ident)[1] %in% 
            unique(object@meta.data[name][, 1])){
            factor <- name
        }
    }
    
    #Marks whether we have to transfer the vect in hidden or not
    hidden <- FALSE
    
    #makes a count table for each cluster
    for (i in 1:length(data_lst)){
        #hidden exists so we mark hidden and skip making a count table for it
        if (i == length(data_lst) & names(data_lst)[i] == "hidden"){
            hidden <- TRUE
        }else {
            #pseudobulks and subsets the cluster
            cnt_tbl <- to_pseudobulk(
                    data_lst[[i]], #The source of what we're generating a count
                    replicate_col = "ratID",
                    cell_type_col = "seurat_clusters",
                    label_col = factor
                    )[[levels(object$seurat_clusters)[match(names(data_lst)[i], 
                       levels(object@active.ident))]]]
            
            #Adds table to list
            if (i != 1){
                furniture <- list.append(furniture, clust_nam = cnt_tbl)
              
            }else {
                #First col = special case as mat is unpopulated
                furniture[[1]] <- cnt_tbl
            }
            
            #Ensures table has correct name
            names(furniture)[i] <- names(data_lst)[i]
        }
    }
    
    if (incl_all){
        #Makes the element with all the data
        allset <- FindClusters(object, resolution = 0)
        Idents(allset) <- factor
        allset <- subset(allset, idents = levels(data_lst[[1]]))
        
        #psuedobulks all and turns into a table to add to list
        all_tbl <- to_pseudobulk(
            allset, #The source of what we're generating a count
            replicate_col = "ratID",
            cell_type_col = "seurat_clusters",
            label_col = factor
        )[["0"]]  #<- The "0" was an inference that works
        furniture <- list.append(furniture, all = all_tbl)
        
        names(furniture)[match("all", names(furniture))] <- paste("all", 
                                  toupper(deparse(substitute(object))), sep = "")
    }
    
    #checks for hidden
    if (hidden){
        #prevents ratID use
        data_lst[["hidden"]] <- data_lst[["hidden"]][data_lst[["hidden"]] != 
                                                                         "ratID"]
        
        #Users can append to hidden. This prevents the code more than 2 factors
        #from existing in hidden when it is used.
        if (length(data_lst[["hidden"]]) > 2){
            data_lst[["hidden"]] <- data_lst[["hidden"]][1:2]
        }
      
        #list that will go in hidden
        hidd_list <- list()
        
        #gets colnames to put in metadata
        for (i in 1:length(data_lst[["hidden"]])){
            #pseudobulks to get colnames
            metas <- colnames(to_pseudobulk(
                data_lst[[1]], #The source of what we're generating a count
                replicate_col = "ratID",
                cell_type_col = "seurat_clusters",
                label_col = data_lst[["hidden"]][i]
                )[[levels(object$seurat_clusters)[
                match(names(data_lst)[1], levels(object@active.ident))]]])
            
            #appends and names 
            hidd_list <- list.append(hidd_list, newCols = metas)
            names(hidd_list)[i] <- data_lst[["hidden"]][i]
        }
        
        #transfers hidden to furniture
        furniture[["hidden"]] <- hidd_list
    }
    
    #returns list of count tables
    return(furniture)
}
