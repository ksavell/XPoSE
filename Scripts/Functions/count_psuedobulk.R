#' Count Psuedobulk
#' 
#' #Makes a list containing the count tables for each of the user's chosen 
#' clusters
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
count_pseudobulk <- function(data_lst, object){
    
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
    
    #makes a count table for each cluster
    for (i in 1:length(data_lst)){
        #pseudobulks and subsets the cluster
        cnt_tbl <- to_pseudobulk(
            data_lst[[i]], #The source of what we're generating a count
            replicate_col = "ratID",
            cell_type_col = "seurat_clusters",
            label_col = factor
        )[[levels(object$seurat_clusters)[i]]]

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
    furniture <- list.append(furniture, "AllClust" = all_tbl)

    #returns list of count tables
    return(furniture)
}
