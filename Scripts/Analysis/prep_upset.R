# I attempt to remake the tibbles with all the genes included

# By Drake Thompson

# Info --------------------------------------------------------------------
#To be filled later

#bin_search function code 



# Loading -----------------------------------------------------------------

library(tibble)


# Function  ---------------------------------------------------------------

#' Prep Upset
#'
#' @param res_lst list of result tibbles
#' @param object Seurat object used to make result tibbles
#'
#' @return tibbles with all genes in our dataset
#' @export
#'
#' @examples
prep_upset <- function(res_lst, object){
    #will be returned at function's end
    return_lst <- list()
    
    ind <- 0
    
    for (cluster in names(res_lst)) {
        ind <- ind + 1
        print(paste("On Index", ind, "of list."))
      
        #allows us to take rows from the base tibble
        cls_frm <- data.frame(res_lst[[cluster]])
        #print(cls_frm)
      
        #makes what will be new tibble and fills genes
        new_tbl <- data.frame(matrix(
                     nrow = length(object@assays[["RNA"]]@counts@Dimnames[[1]]),
                     ncol = ncol(res_lst[[cluster]])))
        
        #Adds some necessary data to new tibble
        new_tbl[, 1] <- object@assays[["RNA"]]@counts@Dimnames[[1]]
        colnames(new_tbl) <- colnames(cls_frm)
        #print(new_tbl)
        
        for (i in 1:nrow(new_tbl)) {
            #print(i)
            if (new_tbl[i, 1] %in% cls_frm[, 1]){
                #replaces row with data in cluster
                new_tbl[i, ] <- cls_frm[match(new_tbl[i, 1], cls_frm[, 1]), ]
            }
        }
        
        #returns to tibble and apppends to list
        new_tbl <- tibble(new_tbl)
        
        return_lst <- list.append(return_lst, new_tbl)
    }
    
    #ensures names are same
    names(return_lst) <- names(res_lst)
    
    #returns new list
    return(return_lst)
}


bin_search <- function(cntnr, target, left = 1, right = length(cntnr)){
  
  if (right >= 2 & left <= right){
    #finds the middle
    #print(right)
    mid <- as.integer((right + left)/2)
    #print(mid)
    
    #found
    if (cntnr[mid] == target){
      #print(mid)
      return(mid)
    
    #target greater
    } else if (cntnr[mid] < target){
      return(bin_search(cntnr, target, mid + 1, right))
      
    #target lesser
    } else {
      return(bin_search(cntnr, target, left, mid - 1))
    }
  }else {
    return(NA)
  }
}


prep_upset_b <- function(res_lst, object){
  #will be returned at function's end
  return_lst <- list()
  
  ind <- 0
  
  #sorted vector of genes
  s_genes <- sort(object@assays[["RNA"]]@counts@Dimnames[[1]], method = "quick")
  
  for (cluster in names(res_lst)) {
    ind <- ind + 1
    print(paste("On Index", ind, "of list."))
    
    #allows us to take rows from the base tibble
    cls_frm <- data.frame(res_lst[[cluster]])
    #print(cls_frm)
    
    #makes what will be new tibble and fills genes
    new_tbl <- data.frame(matrix(
      nrow = length(object@assays[["RNA"]]@counts@Dimnames[[1]]),
      ncol = ncol(res_lst[[cluster]])))
    
    #Adds some necessary data to new tibble
    new_tbl[, 1] <- object@assays[["RNA"]]@counts@Dimnames[[1]]
    colnames(new_tbl) <- colnames(cls_frm)
    rownames(new_tbl) <- new_tbl[, "gene"]
    #print(new_tbl)
    
    for (i in 1:nrow(new_tbl)) {
      #print(i)
      if (!is.na(bin_search(cls_frm[, 1], new_tbl[i, 1]))){
      #if (new_tbl[i, 1] %in% cls_frm[, 1]){
        #replaces row with data in cluster
        #print(bin_search(s_genes, new_tbl[i, 1]))
        new_tbl[i, ] <- cls_frm[s_genes[bin_search(s_genes, new_tbl[i, 1])], ]
      }
    }
    
    #returns to tibble and apppends to list
    new_tbl <- tibble(new_tbl)
    print(new_tbl)
    
    return_lst <- list.append(return_lst, new_tbl)
  }
  
  #ensures names are same
  names(return_lst) <- names(res_lst)
  
  #returns new list
  return(return_lst)
}




