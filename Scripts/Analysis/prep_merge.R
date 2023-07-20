# I attempt to remake the tibbles with all the genes included

# By Drake Thompson

# Info --------------------------------------------------------------------
# This contains steps to merge the results tables for things like upset plots. 
# This utilizes a binary search algorithm which means that the data needs to be 
# sorted. A version of the function that doesn't require the data to be sorted 
# is included here. It runs ~3-7 seconds slower per cluster.

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
    
    #goes through entire results list
    for (cluster in names(res_lst)) {
        #Tells user where function is in correction process
        ind <- ind + 1
        print(paste("On Index ", ind, "/", length(res_lst), sep = ""), quote = FALSE)
        
        #allows us to take rows from the base tibble
        cls_frm <- data.frame(res_lst[[cluster]])
        
        #new tibble that will added to list
        new_tbl <- data.frame(matrix(
                     nrow = length(object@assays[["RNA"]]@counts@Dimnames[[1]]),
                     ncol = ncol(res_lst[[cluster]])))
        
        #Adds necessary data to new tibble
        new_tbl[, 1] <- object@assays[["RNA"]]@counts@Dimnames[[1]]
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




#' Merge Results
#' 
#' This merges result tables generated by DEseq. It can further filter the data
#' based on a user's chosen significance number. 
#'
#' @param agg_lst a list of results tibbles generated by DEseq
#' @param filter TRUE/FALSE that tells function whether you want to filter out
#'               certain rows. Set to true by default
#' @param thres a value to check for significance with. Default is 0.05
#'
#' @return a merged result table
#' @export
#'
#' @examples
merge_results <- function(agg_lst, filter = TRUE, thres = 0.05){
    #checks length
    if (length(agg_lst) > 1){
        
        #error marker
        diff_cnts_flag <- FALSE
        
        #checks to see is all row counts are the same
        for (clustA in agg_lst[1:(length(agg_lst) - 1)]) {
            for (clustB in agg_lst[2:length(agg_lst)]){
                if (nrow(clustA) != nrow(clustB)){
                  diff_cnts_flag <- TRUE
                }
            }
        }
        
        #error thrown to stop potential errors
        if (diff_cnts_flag){
            stop("List '", deparse(substitute(agg_lst)), "' has tibbles with ",
                 "different row counts.\n  ", 
                 "Please use a function like 'prep_merge' to ensure rows are all same",
                 "number\n")
        }
        
        #modifies list 
        ag <- agg_lst
        
        #list performs a bit funky with odd numbered lists, this prevents small issues
        #that can occur because of that
        odd <- FALSE
        if (length(ag) %% 2 == 1){
            ag[["dupe"]] <- ag[[1]]
            odd <- TRUE
        }
        
        #loops through list and merges
        for (i in 1:(length(ag) - 1)){
            ag[[2]] <- merge(ag[[1]], ag[[2]], 
                             by='gene', all=TRUE, 
                             no.dups = TRUE,
                             suffixes=c(paste("_", names(ag)[1], sep = ""),
                                        paste("_", names(ag)[2], sep = "")))
            ag <- ag[-1]
        }
        
        #removes extra column created so merging isn't funky
        if(odd){
            ag[[1]] <- ag[[1]][!(str_detect(names(ag[[1]]), "dupe"))]
        }
        
        aggr <- ag[[1]]
    }
    else {
        aggr <- as.data.frame(agg_lst[[1]])
    }
  
    #marks rows 
    rownames(aggr) <- ag[[1]][, "gene"]
    aggr <- aggr[, colnames(aggr) != "gene"] 
    
    #filters data if user wants it to be done
    if (filter){
        #creates new column
        aggr$keep <- FALSE
        
        #marks significant rows as TRUE
        for (name in names(aggr)[str_detect(names(aggr), "padj")]){
            aggr$keep[aggr[, name] < thres] <- TRUE
        }
        
        #removes unmarked rows
        aggr <- aggr[aggr$keep, ]
        
        #small check
        cat("There should be a TRUE under this line\n",
            !(FALSE %in% aggr$keep), "\n", sep = "")
        
        #removes extra col
        aggr <- aggr[, colnames(aggr) != "keep"] 
    }
    
    #returns the big table as data frame
    return(aggr)
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

