#make an upset without user input




#All functions

#' Make Upset
#' 
#' Desc
#'
#' @param seur_obj1 
#' @param seur_obj2 
#' @param factor 
#' @param comp_vect 
#' @param inc_all 
#'
#' @return upset plot data
#' @export
#'
#' @examples
make_upset <- function(seur_obj1, seur_obj2, threshold, factor, comp_vect, incl_all = FALSE, genes_sorted = FALSE){
  #libraries
  library(dplyr)
  library(purrr)
  library(Seurat)
  library(Libra)
  library(here)
  library(rlist)
  library(stringr)
  library(tibble)
  library(DESeq2)
  
  #a series of error messages
  #code last
  
  #list of objects for abstraction purposes
  obj_list <- list(obj1 = seur_obj1, obj2 = seur_obj2)
  names(obj_list) <- c(deparse(substitute(seur_obj1)), 
                       deparse(substitute(seur_obj1)))
  
  #list of merged objects to be merged into co-expression table
  merg_lst <- list()
  
  for (obj in obj_list){
    #prep pseudobulk
    #prep_pseudobulk(obj, threshold, factor, comp_vect)
    
    #print(obj)
    #count pseudobulk tbl lst
    
    
    #prep DESeq
    
    #run DEseq
    execute_DEseq(run_pseudobulk(obj, threshold, factor, comp_vect, incl_all), 
                  obj, factor, comp_vect)
    #prep merge
    
    #merge
    
    #append to merg_lst
  }
  
  #make the inter and distinct mats for all in merg
  
  #return lst that has 4 elems (disUp, disDown, interUp, interDown)
  
}


#' Verify Factor
#' 
#' General function that will take in a vector of user input and correct it using
#' a vector of the wanted factor's names. 
#'
#' @param object SE
#' @param vect the vector being verified, should have at least 2 elems
#' @param factor the factor to verify with. Within our object this will be 
#'               'group' or 'sex'
#'
#' @return the user's (corrected) choice vector
#' @export
verify_factor <- function(object, vect, factor){
  #makes vector of valid options, normal and lowercase
  valids <- unique(object[[factor]])[, 1]
  valids_low <- tolower(valids)

  #There is a case where factor is the first elem so this ensures output to user
  #is consistent
  mod <- 0

  #Loops vector
  for (i in 1:length(vect)){
    #Checks if first elem is the factor
    if (i == 1){
      if (vect[i] != factor){
        #Loops until valid input
        while (!(vect[i] %in% valids) & vect[i] != "d"){
          #checks for case, converts to proper if not
          if (tolower(vect[i]) %in% valids_low){
            vect[i] <- valids[match(tolower(vect[i]), valids_low)]
          }else {
            #prompts user about
            cat("Invalid input for Comparison", i,
                ". Please check for typo.\n",
                "Valid inputs include:\n",
                "If you would like to delete this option, please enter 'd'. \n",
                sep = "")
            sapply(valids, print, quote = FALSE)
            vect[i] = readline(paste("Comparsion", i,
                                     ": ", sep = ""))
          }
        }
        #If the factor is the first index, i needs to be 1
      } else {
        mod <- 1
      }
    }else {
      #Loops until valid input
      while (!(vect[i] %in% valids) & vect[i] != "d"){
        #checks for case
        if (tolower(vect[i]) %in% valids_low){
          vect[i] <- valids[match(tolower(vect[i]), valids_low)]
        }else {
          cat("Invalid input for Comparison", i - mod,
              ". Please check for typo.\n",
              "Valid inputs include:\n",
              "If you would like to delete this option, please enter 'd'. \n",
              sep = "")
          sapply(valids, print, quote = FALSE)
          vect[i] = readline(paste("Comparsion", i - mod, ": ", sep = ""))
        }
      }
    }
  }

  #Check to see if we need to rerun the function
  rerun <- FALSE

  #dupe check
  if (length(vect) > 1){
    for (i in (1 + mod):(length(vect) - 1)){
      for(j in (2 + mod):(length(vect))){
        #Prevents issue where i and j can both be 2
        if (i != j){
          while (vect[i] == vect[j]){
            rerun <- TRUE
            cat("Dupe found at Comparsion", j - mod,
                ". Please use a different option within your set.\n",
                "Valid inputs include:\n",
                "If you would like to delete this option, please enter 'd'. \n",
                sep = "")
            sapply(valids, print, quote = FALSE)
            cat("Please don't put in the same thing in...")
            vect[i + 1] = readline(paste("Comparsion", j - mod, ": ", sep = ""))
          }
        }
      }
    }
    #Cleans vect of d's
    vect <- vect[vect != "d"]
  }

  while (length(vect) < 2){
    rerun <- TRUE
    cat("Additional comparison needed.",
        " Please enter one of the below options:\n")
    sapply(valids, print, quote = FALSE)
    vect <- append(vect, readline("New Comparsion: "))
  }

  #Validates possible new changes
  if (rerun){
    vect <- verify_factor(object, vect, factor)
  }

  #returns vect with changes (if any)
  return(vect)
}


#' Run Psuedobulk
#' 
#' gets clusters and creates count tables for them
#'
#' @param object 
#' @param threshold 
#' @param factor 
#' @param comp_vect 
#' @param incl_all 
#'
#' @return
#' @export
#'
#' @examples
run_pseudobulk <- function(object, threshold, factor, comp_vect, incl_all){
  #Makes base data table so we can import the data to the frame
  data_tbl <- table(object$seurat_clusters, object[[factor]][, 1])
  
  #Makes data frame to store
  clust_tab <- data.frame(matrix(ncol = length(unique(object[[factor]][, 1])), 
                                 nrow = length(rownames(data_tbl))))
  
  #labels the frame's rows and columns
  colnames(clust_tab) <- unique(object[[factor]][, 1])
  rownames(clust_tab) <- unique(levels(object@active.ident))
  
  #populates data frame
  for (i in 1:length(unique(object[[factor]][, 1]))){
    clust_tab[, unique(object[[factor]][, 1])[i]] <- data_tbl[, i]
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
  #Makes cluster subsets derived from user input
  for (i in 1:length(incl_clust)) {
    
    #Makes val of data_lst the subset
    data_lst[[incl_clust[i]]] <- subset(object, idents = incl_clust[i])
    
    #Splitting stuff
    Idents(data_lst[[incl_clust[i]]]) <- factor
    data_lst[[incl_clust[i]]] <- subset(data_lst[[incl_clust[i]]], 
                                   idents = comp_vect)
  }
 
  #returns list
  #return(data_lst)
  #print(data_lst)
  
  furniture <- list(matrix(1:length(object@assays[["RNA"]]@counts@Dimnames[[1]]),
                           ncol = 1))
  
  #makes a count table for each cluster
  for (i in 1:length(data_lst)){
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


execute_DEseq <- function(tbl_lst, object, factor, comp_vect){
  
  #make metadata, call it met_fram
  #splits names into list of the ratID's and groups
  temp <- sapply(names(tbl_lst[[1]]), strsplit, ":")
  justF <- c()
  
  #Fills justF w/ only the comparisons in the split
  for (i in 1:length(temp)){
    justF <- append(justF, temp[[i]][2]) 
  }
  
  #Makes data frame and labels it 
  met_fram <- data.frame(justF)
  rownames(met_fram) <- names(tbl_lst[[1]])
  colnames(met_fram) <- factor 
  
  #checks that the metadata is set up correctly
  cat("There should be at least one 'TRUE' under this line:\n")
  cat(as.character(unique(met_fram[,1]) %in% unique(object@meta.data[factor][, 1])), "\n")
  cat("\n")
  
  #returns frame
  #return(met_fram)
  
  #makes necessary materials from params
  contrast_vect <- c(factor, comp_vect)
  form_fact <- as.formula(paste("~", factor, sep = ""))
  
  #our base list to hold tibbles, will be returned at function's end
  agg_res <- list()
  
  #Makes result table for each count table in list (excludes hidden)
  for (i in 1:length(tbl_lst)){
      #Makes object
      cluster <- DESeqDataSetFromMatrix(tbl_lst[[i]],  #Data Table
                                        colData = met_fram,  #metaData
                                        design = form_fact)
      
      #Filters out insignificant parts
      cluster <- cluster[ rowSums(counts(cluster)) > 5, ]
      cluster <- DESeq(cluster)
      
      #Makes results obj
      results <- results(cluster,
                         contrast = contrast_vect,
                         alpha = 0.05,
                         pAdjustMethod = "fdr")
      
      #Turns results into a Tibble
      results_tib <- results %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
      
      #appends to list of tibbles
      agg_res <- list.append(agg_res, clust_nam = results_tib)
      names(agg_res)[i] <- names(tbl_lst)[i]
  }
  #returns results
  print(agg_res)
  return(agg_res)
}






#make_upset(glut, gaba, 100,  "group", c("positive", "negative"))
