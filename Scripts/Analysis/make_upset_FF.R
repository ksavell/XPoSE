#make an upset without user input




#All functions

#' Make Upset
#' 
#' Will return 4 tables to the user with inputted parameters so that an Upset plot
#' can be easily made in Prism.
#'
#' @param seur_obj1 one of the Seurot objects to do analysis on. Order doesn't 
#'                  matter
#' @param seur_obj2 Same as above
#' @param threshold minimum number of nuclei per comparison to include a cluster 
#'                  in analysis
#' @param factor An aspect of your data like population or sex that you will use 
#'               to anaylze your data
#' @param comp_vect Aspects within the factor you would like to compare. First 
#'                  elem will become your Experimental Variable, the second will
#'                  become your control
#' @param p_sig Significance Value to filter by when deciding which rows to keep
#'              when merging. Default is 0.05. 
#' @param incl_all logical, determines whether you want to include an all 
#'                 clusters variable in analysis. Default is FALSE.
#' @param genes_sorted logical, marks if the faster search algorithm is usable.
#'                     If marked as true and the gene list still isn't sorted, 
#'                     it will act as if this was marked false. Default is FALSE.
#'
#' @return upset plot data
#' @export
#'
#' @examples
make_upset <- function(seur_obj1, seur_obj2, threshold, factor, comp_vect, 
                       p_sig = 0.05, incl_all = FALSE, genes_sorted = FALSE){
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
  library(UpSetR)
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(data.table)
  
  #for process time checking
  timer <- proc.time()
  
  #a series of error messages
  #code last
  
  #list of objects for abstraction purposes
  obj_list <- list(obj1 = seur_obj1, obj2 = seur_obj2)
  names(obj_list) <- c(deparse(substitute(seur_obj1)), 
                       deparse(substitute(seur_obj1)))
  
  #In the case that the gene lists are not the same, this ensures the merge 
  #doesn't cause any errors
  gene_vect <- c()
  gene_lngth_s <- TRUE
  print("eeeee")
  #checks for different lengths or different genes within
  if ((length(seur_obj1@assays[["RNA"]]@counts@Dimnames[[1]]) != 
       length(seur_obj2@assays[["RNA"]]@counts@Dimnames[[1]])) |
      (FALSE %in% (seur_obj1@assays[["RNA"]]@counts@Dimnames[[1]] == 
                  seur_obj2@assays[["RNA"]]@counts@Dimnames[[1]]))){
    print("wa")
    #gene vectors are not same length
    gene_lngth_s <- FALSE
    
    #make the unique of c(obj1, obj2) and sorts 
    gene_vect <- sort(unique(c(seur_obj1@assays[["RNA"]]@counts@Dimnames[[1]], 
                               seur_obj2@assays[["RNA"]]@counts@Dimnames[[1]])), 
                      method = "quick")
    
    #return(gene_vect)
  }
  
  #list of merged objects to be merged into co-expression table
  merg_lst <- list()
  
  ind <- 1
  #pseudobulks, runs DESeq and merges tables
  for (obj in obj_list){
    #differing lengths or genes within
    if (gene_lngth_s){
      gene_vect <- obj@assays[["RNA"]]@counts@Dimnames[[1]]
    }
    
    #gets results tables
    obj_rt <- execute_DEseq(run_pseudobulk(obj, threshold, factor, comp_vect, incl_all), 
                  obj, factor, comp_vect)
    
    #checks if data sorted if genes_sorted is TRUE
    if (genes_sorted & 
        !(FALSE %in% (obj@assays[["RNA"]]@counts@Dimnames[[1]] == 
                      sort(obj@assays[["RNA"]]@counts@Dimnames[[1]], 
                           method = "quick")))){
      #runs version w/ binary search and merges
      #obj_rt <- merge_results(prep_merge_fast(obj_rt, obj, gene_vect), TRUE, p_sig)
      obj_rt <- prep_merge_fast(obj_rt, obj, gene_vect)
      
    #List unsorted or user did not specify
    }else {
      #run unsorted version and merges
      obj_rt <- prep_merge(obj_rt, obj, gene_vect)
    }
    
    #return(obj_rt)
    #append to merg_lst and name
    merg_lst <- list.append(merg_lst, obj_rt)
    names(merg_lst)[ind] <- paste("thing", ind, sep = "")
    
    ind <- ind + 1
    #return(merg_lst)
  }
  
  merg_lst <- append(merg_lst[[1]], merg_lst[[2]])
  #return(merg_lst)
  
  coexp <- merge_results(merg_lst, thres = p_sig)
  #return(coexp)
  #return(merg_lst)
  #merge the objs
  
  #make up and down data frames
  
  ##Paralelel vectors are your friend
  
  #special str detect for loop, needs to store names in another vector for later 
  #use
  
  print(names(coexp))
  
  #make the inter and distinct mats for all in merg
  
  #return lst that has 4 elems (disUp, disDown, interUp, interDown)
  
  #tells user how long function took
  elpd <- (proc.time() - timer)[[3]]
  print(paste("Process took", round(elpd, digits = 3), "seconds"), 
        quote = FALSE)
  
  return(coexp)
}


#' Verify Factor
#' 
#' General function that will take in a vector of user input and correct it using
#' a vector of the wanted factor's names. 
#'
#' @param object the Seurot object to do analysis on
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


#' Title
#'
#' @param tbl_lst a list of count tables generated from cluster data
#' @param object the Seurot object to do analysis on
#' @param factor the factor to verify with. Within our object this will be 
#'               'group' or 'sex' 
#' @param comp_vect Aspects within the factor you would like to compare. First 
#'                  elem will become your Experimental Variable, the second will
#'                  become your control
#'
#' @return a list of results tables generated from the count tables
#' @export
#'
#' @examples
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
  #print(agg_res)
  return(agg_res)
}


#' Prep Merge
#' 
#' Ensures all tibbles in result list align so they can be merged. Doesn't 
#' requires sorted data so it is a bit slower.
#'
#' @param res_lst 
#' @param object 
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
    new_tbl <- data.frame(matrix(
      nrow = length(gene_vect),
      ncol = ncol(res_lst[[cluster]])))
    
    #Adds some necessary data to new tibble
    new_tbl[, 1] <- gene_vect
    colnames(new_tbl) <- colnames(cls_frm)
    rownames(cls_frm) <- cls_frm[, "gene"]
    
    #print("Vgf" %in% gene_vect)
    #print(match("Vgf",  new_tbl[, 1]))
    #print(new_tbl[28708 , 1] == "Vgf")
    #print(28708 %in% which(new_tbl[, 1] %in% res_lst[[cluster]]$gene))
    #print(which(res_lst[[cluster]]$gene %in% new_tbl[, 1]))
    #for (i in which(res_lst[[cluster]]$gene %in% new_tbl[, 1])) {
    #for (i in 1:length(res_lst[[cluster]]$gene %in% object@assays[["RNA"]]@counts@Dimnames[[1]])) {
    for (i in which(new_tbl[, 1] %in% res_lst[[cluster]]$gene)) {
      # if (new_tbl[i , 1] == "Vgf"){
      #   print(cls_frm[new_tbl[i, "gene"], ])
      #   print(new_tbl[i, "gene"])
      # }
      #print(new_tbl[i, "gene"])
      #replaces row with data in cluster
      #new_tbl[i, ] <- cls_frm[new_tbl[i, "gene"], ]
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


#' Prep Merge (Fast)
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
prep_merge_fast <- function(res_lst, object, gene_vect = object@assays[["RNA"]]@counts@Dimnames[[1]]){
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
      nrow = length(gene_vect),
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
  #modifies list 
  ag <- agg_lst
  
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
  rownames(aggr) <- as.data.frame(ag[[1]])[, "gene"]
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


#make_upset(glut, gaba, 100,  "group", c("positive", "negative"))
