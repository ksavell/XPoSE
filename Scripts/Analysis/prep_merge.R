# I attempt to remake the tibbles with all the genes included

# By Drake Thompson

# Info --------------------------------------------------------------------
#To be filled later

#bi_search function code modeled after the recursive search function by Geeks
#for geeks. Link here:
#https://www.geeksforgeeks.org/binary-search/



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
prep_upset_unsrtd <- function(res_lst, object){
    timer <- proc.time()
    #will be returned at function's end
    return_lst <- list()
    
    ind <- 0
    
    for (cluster in names(res_lst)) {
        ind <- ind + 1
        print(paste("On Index", ind, "of list."))
      
        #allows us to take rows from the base tibble
        cls_frm <- data.frame(res_lst[[cluster]])
      
        #makes what will be new tibble and fills genes
        new_tbl <- data.frame(matrix(
                     nrow = length(object@assays[["RNA"]]@counts@Dimnames[[1]]),
                     ncol = ncol(res_lst[[cluster]])))
        
        #Adds some necessary data to new tibble
        new_tbl[, 1] <- object@assays[["RNA"]]@counts@Dimnames[[1]]
        colnames(new_tbl) <- colnames(cls_frm)
        
        for (i in 1:nrow(new_tbl)) {
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
    
    print(proc.time() - timer)
    
    #returns new list
    return(return_lst)
}

# user  system elapsed 
# 46.499  29.062  75.912 


#' Binary Search
#' 
#' Recursively finds target elem's index within the container using a binary 
#' search algorithm
#' If we release this in a package, I will need to verify it's an atomic vector
#'
#' @param cntnr an atomic vector containing the target
#' @param target elem of cntnr you're trying to find index of
#' @param left first index of wanted range (default is 1)
#' @param right last index of wanted range (default is length(cntnr))
#'
#' @return the index of the target
#' @export
#'
#' @examples
#' #Find index of "O" in LETTERS
#' bi_search(LETTERS, "O")
#' #Find index of "O" in a specific range (27-53)
#' bi_search(c(LETTERS, LETTERS, LETTERS), "O", 27, 53)
bi_search <- function(cntnr, target, left = 1, right = length(cntnr)){
  
  if (right >= 2 & left <= right){
    #finds middle
    mid <- as.integer((right + left)/2)
    
    #target found
    if (cntnr[mid] == target){
      return(mid)
    
    #target greater
    } else if (cntnr[mid] < target){
      return(bi_search(cntnr, target, mid + 1, right))
      
    #target lesser
    } else {
      return(bi_search(cntnr, target, left, mid - 1))
    }
    
  #object not found
  }else {
    return(NA)
  }
}


#' Prep Merge
#' 
#' Ensures all tibbles in result list align so they can be merged. Requires 
#' sorted data as it uses a binary search algorithm.
#'
#' @param res_lst a list of result tibbles
#' @param object Seurat object analysis was performed on
#'
#' @return the corrected result list
#' @export
#'
#' @examples
#' #Prepares tibbles in res object.
#' res <- prep_merge(res, gs)
prep_merge <- function(res_lst, object){
    timer <- proc.time()
    #will be returned at function's end
    return_lst <- list()
    
    ind <- 0
    
    for (cluster in names(res_lst)) {
        ind <- ind + 1
        print(paste("On Index ", ind, "/", length(res_lst), sep = ""), quote = FALSE)
        
        #allows us to take rows from the base tibble
        cls_frm <- data.frame(res_lst[[cluster]])
        
        #makes what will be new tibble and fills genes
        new_tbl <- data.frame(matrix(
                      nrow = length(object@assays[["RNA"]]@counts@Dimnames[[1]]),
                      ncol = ncol(res_lst[[cluster]])))
        
        #Adds some necessary data to new tibble
        new_tbl[, 1] <- object@assays[["RNA"]]@counts@Dimnames[[1]]
        colnames(new_tbl) <- colnames(cls_frm)
        rownames(cls_frm) <- cls_frm[, "gene"]
        
        #fills shared rows w/ base data
        for (i in 1:nrow(new_tbl)) {
            if(!is.na(bi_search(cls_frm[, 1], new_tbl[i, 1]))){
                #replaces row with data in cluster
                new_tbl[i, ] <- cls_frm[bi_search(cls_frm[, 1], new_tbl[i, 1]), ]
            }
        }
        
        #returns to tibble and apppends to list
        new_tbl <- tibble(new_tbl)
        
        return_lst <- list.append(return_lst, new_tbl)
    }
    
    #ensures names are same
    names(return_lst) <- names(res_lst)
    
    el <- (proc.time() - timer)[[3]]
    print(paste("Process took", round(el, digits = 3), "seconds"), quote = FALSE)
    
    #returns new list
    return(return_lst)
}


# user  system elapsed 
# 37.593  22.420  60.311 


prep_upset_c <- function(res_lst, object){
  timer <- proc.time()
  #will be returned at function's end
  return_lst <- list()
  
  ind <- 0
  
  for (cluster in names(res_lst)) {
    ind <- ind + 1
    print(paste("On Index", ind, "of list."))
    
    #allows us to take rows from the base tibble
    cls_frm <- data.frame(res_lst[[cluster]])
    
    #makes what will be new tibble and fills genes
    new_tbl <- data.frame(matrix(
      nrow = length(object@assays[["RNA"]]@counts@Dimnames[[1]]),
      ncol = ncol(res_lst[[cluster]])))
    
    #Adds some necessary data to new tibble
    new_tbl[, 1] <- object@assays[["RNA"]]@counts@Dimnames[[1]]
    colnames(new_tbl) <- colnames(cls_frm)
    #rownames(cls_frm) <- cls_frm[, "gene"]
    #rownames(new_tbl) <- new_tbl[, "gene"]
    
    for (i in 1:length(cls_frm[, 1])){
      new_tbl[bi_search(new_tbl[, 1], cls_frm[i, 1]), ] <- cls_frm[i, ]
    }
    
    #returns to tibble and apppends to list
    new_tbl <- tibble(new_tbl)
    
    return_lst <- list.append(return_lst, new_tbl)
  }
  
  #ensures names are same
  names(return_lst) <- names(res_lst)
  
  print(proc.time() - timer)
  
  #returns new list
  return(return_lst)
}




# timer <- proc.time()
# for (gene in bbb){
# match(gene, bbb)
# }
# proc.time() - timer
# 
# cat("\n\n\n\n\n\n\n\n")
# 
# timer <- proc.time()
# for (gene in bbb){
# bi_search(bbb, gene)
# }
# g <- (proc.time() - timer)[[3]]

