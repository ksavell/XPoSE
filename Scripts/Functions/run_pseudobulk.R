#' Run Psuedobulk
#' 
#' gets clusters and creates count tables for them
#'
#' @param object the Seurot object to do analysis on
#' @param threshold minimum number of nuclei per comparison to include a cluster 
#'                  in analysis
#' @param factor An aspect of your data like population or sex that you will use 
#'               to anaylze your data
#' @param comp_vect Aspects within the factor you would like to compare. First 
#'                  elem will become your Experimental Variable, the second will
#'                  become your control
#' @param incl_all logical, determines whether you want to include an all 
#'                 clusters variable in analysis
#'
#' @return count tables for each chosen cluster
#' @export
#'
#' @examples
run_pseudobulk <- function(object, threshold, factor, comp_vect, incl_all){
    #Makes base data table so we can import the data to the frame
    data_tbl <- table(object$cluster_name, object[[factor]][, 1])
    
    #Makes data frame to store
    clust_tab <- data.frame(matrix(ncol = length(unique(object[[factor]][, 1])), 
                                   nrow = length(rownames(data_tbl))))
    
    #labels the frame's rows and columns
    colnames(clust_tab) <- colnames(data_tbl)
    rownames(clust_tab) <- rownames(data_tbl)
    
    #populates data frame
    for (i in 1:length(colnames(clust_tab))){
      clust_tab[, colnames(clust_tab)[i]] <- data_tbl[, i]
    }
    
    #Prints table for user pleasure
    cat("\nTable:\n")
    final <- clust_tab[clust_tab[, comp_vect[1]] > threshold &
                         clust_tab[, comp_vect[2]] > threshold, ] %>% 
                         select(comp_vect[1], comp_vect[2])
    
    print(final)
    cat("\n")
    
    incl_clust <- rownames(final)
    
    data_lst <- list()
    Idents(object) <- "cluster_name"
    #Makes cluster subsets derived from user input
    for (i in 1:length(incl_clust)) {
        
        #Makes val of data_lst the subset
        data_lst[[rownames(final)[i]]] <- subset(object, idents = rownames(final)[i])
        
        #Splitting stuff
        Idents(data_lst[[incl_clust[i]]]) <- factor
        data_lst[[incl_clust[i]]] <- subset(data_lst[[incl_clust[i]]], 
                                       idents = comp_vect)
    }
    
    furniture <- list(matrix(1:length(object@assays[["RNA"]]@counts@Dimnames[[1]]),
                             ncol = 1))
    
    #makes a count table for each cluster
    for (i in 1:length(data_lst)){
        #pseudobulks and subsets the cluster
        cnt_tbl <- to_pseudobulk(
          data_lst[[i]], #The source of what we're generating a count
          replicate_col = "ratID",
          cell_type_col = "cluster_name",
          label_col = factor
        )[[unique(object$cluster_name)[match(names(data_lst)[i], 
                                                unique(object$cluster_name))]]]
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
    
    #print(furniture)
    return(furniture)
}
