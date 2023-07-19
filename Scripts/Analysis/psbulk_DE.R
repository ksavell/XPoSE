# Remakes the remade base file with code sections, minor improvements and such

# By: Drake Thompson, Padmashri Saravanan, Katherine Savell

# Info --------------------------------------------------------------------

# This script creates: 
#           ** a 'MEGA' table that can divided into subtables 
#                       (use prep_DEseq to get subtables)
#           ** metadata table (stored  as metaData) for making result tables
#                       (use get_result_tbl or run_DEseq)

# Making individual count and result tables will have to be done by user
# This is done with the three remaining functions:
#                                              prep_DEseq()
#                                              get_result_tbl()
#                                              run_DEseq()

# Note: This code is broken into various sections. You can collapse sections and 
# look at specific chunks of code. There is also a Jump To Command that takes you 
# to the specific sections under 'Code' on the top menu bar.


# Loading -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
# loads all required packages
library(dplyr)
library(purrr)
library(Seurat)
library(Libra)
library(here)


## Data Loading ------------------------------------------------------------
# Loads data depending on user choice
prompt <- paste("Use both data sets?\n", 
                "1. Yes\n",
                "2. Just Glut\n",
                "3. Just GABA\n")
cat(prompt)
choice = readline()

# validates input, loops if invalid
while (!(choice %in% 1:3)) {
  # "error" message
  print("Invalid Input. Please return a 1, 2, or 3")
  
  # asks fo input again
  cat(prompt)
  choice = readline()
}

# Loads data, you may need to change below code depending on device. 
##-----------------------------------------------##
# Loads data depending on user choice
if (choice == 1 | choice == 2){
  to_load = readline("Load glut? (y/n): ")
  if (to_load == "y"){
    load("data/glut_subset_forthpipe_11102022.RData") # change to here() 
    glut <- glut.subset2
  }
}
if (choice == 1 | choice == 3){
  to_load = readline("Load GABA? (y/n): ")
  if (to_load == "y"){
    #load("data/gaba_subset_forthpipe_03142023.RData") #Newer doesn't work with this :P
    load("data/gaba_subset_forthpipe_11102022.RData")
    gaba <- gaba.subset2
  }
}
##------------------------------------------------##


# Function Section # --------------------------------------------------------
## Notes About Function Section --------------------------------------------
# SE refers to 'self explanatory'
# Functions are collapsable if you want to just read descriptions and this section
# feels cluttered.

## List of Functions that May Require Changes ------------------------------
##------------------------##
# List of functions that may need to be changed with a new experiment
##------------------------##
### verify_factor:
#               This takes from a rather specific part of the data, where it pulls
#               from may have to change if data is formatted differently
#
###filter_clusters:
#                  3rd Line: if you want to compare more than 2 groups, you will
#                            will have to add more readlines
#                  20th Line: This will have to be changed if there's a greater
#                             group size
#
###make_cluster_list:
#                  Lines 4/5: If you want more than 2 groups, changes will be 
#                             needed. Things like changing readline and prompt
#
###get_cnt_tbls:
#                  Line 3: If the data is stored differently or the path is diff-
#                          erently, you may need to change this line. It is 
#                          getting the length of the amount of gene names. 
#                  Line 13: If you want more than the groups (say you also want 
#                           sex), you may need to add that. But I'm unsure. 
#
###make_meta: 
#           The whole function may need to be reworked a bit if you want to store 
#           more than just 'group' in metadata. Say if you also needed 'sex' for 
#           the metadata. Otherwise should work whenever you just need 'group'.
#
###prep_DEseq: 
#              The basis of this function findind the piece is on ratID so if its
#              location changes in the object, problems may occur. Please be aware
#              of that!





##Functions -----------------------------------


#' Verify Factor
#' 
#' General function that will take in a vector of user input and correct it using
#' a vector of the wanted factor's names. 
#'
#' @param object SE
#' @param vect the vector being verified
#' @param factor the factor to verify with. Within our object this will be 
#'               'group' or 'sex'
#'
#' @return the user's (corrected) choice vector
#' @export
verify_factor <- function(object, vect, factor){
  
  #makes vector of valid options
  valids <- unique(object[[factor]])[,1]
  
  #Loops vector
  for (i in 1:length(vect)){
    #There is a case where factor is the first elem so this ensure output to user
    #is consistent
    mod <- 0 
    if (i == 1){
      if (vect[i] != factor){
        #Loops until valid input
        while (!(vect[i] %in% valids)){
          cat("Invalid input for Comparison", i, ". Please check for typo.\n",
              "Valid inputs include:\n", sep = "")
          sapply(valids, print, quote = FALSE)
          vect[i] = readline(paste("Comparsion", i, ": ", sep = ""))
        }
        #If the factor is the first index, i needs to be 1
      } else {
        mod <- 1
      }
    }else {
      #Loops until valid input
      while (!(vect[i] %in% valids)){
        cat("Invalid input for comparison", (i + mod), ". Please check for typo.\n",
            "Valid inputs include:\n", sep = "")
        sapply(valids, print, quote = FALSE)
        vect[i] = readline(paste("comparison", (i + mod), ": ", sep = ""))
      }
    }
  }
  
  #returns vect with changes (if any)
  return(vect)
}


#' Get Count Tables
#'
#' @param data_lst list of searot objects for dataser, generated by 
#'                 'make_cluster_list'
#' @param clust_vect vector of clusters for object, generated by 'Kowalski'
#' @param object SE, used to make allset
#' @param dataname name of data, used to make all_tbl
#'
#' @return the "MEGA" table that contains all the count tables
#' @export
#'
#' @examples
#' #Makes count tables for the requested clusters for glut object
#' mega_tbl <- get_cnt_tbls(glut_lst, glut_dat, glut, "Glut")
get_cnt_tbls <- function(data_lst, clust_vect, object, dataname){
  
  #Will hold all the tables
  furniture <- matrix(1:length(data_lst[[1]]@assays[["RNA"]]@counts@Dimnames[[1]]), ncol = 1)
  
  #makes a count table for each cluster
  for (i in 1:length(data_lst)){
    #pseudobulks and temporarily stores
    bulk_var <- to_pseudobulk(
      data_lst[[i]], #The source of what we're generating a count
      replicate_col = "ratID",
      cell_type_col = "seurat_clusters", #this is required, could merge as 1 to run
      label_col = "group"
    )
    
    #subset the cluster to test, here they are combined into one called "0"
    cnt_tbl <- bulk_var[[as.character(clust_vect[i])]]
    
    #Makes 'dummy' col and names to mark it
    if (i != 1){
      furniture <- cbind(furniture, 1:nrow(furniture))
      names(furniture)[ncol(furniture)] <- as.character(names(data_lst)[i])
      #First col = special case as mat is unpopulated
    }else {
      colnames(furniture)[1] <- as.character(names(data_lst)[i])
    }
    
    #adds new data to furniture
    furniture <- cbind(furniture, cnt_tbl)
  }
  
  #Makes the element with all the data
  allset <- FindClusters(object, resolution = 0)
  Idents(allset) <- 'group'
  allset <- subset(allset, idents = levels(data_lst[[1]])) #This may not be replicable but it works here
  
  #pseudobulks the all 
  all_dat <- to_pseudobulk(
    allset, #The source of what we're generating a count
    replicate_col = "ratID",
    cell_type_col = "seurat_clusters", #this is required, could merge as 1 to run
    label_col = "group"
  )
  
  #Turns that all into a table
  all_tbl <- all_dat[["0"]]  #<- The "0" was an inference so it may not work elsewhere
  
  #Names sect and adds to mega table
  furniture <- cbind(furniture, 1:nrow(furniture))
  names(furniture)[ncol(furniture)] <- paste("all", dataname, sep = "")
  furniture <- cbind(furniture, all_tbl)
  
  #returns Mega Table, the IKEA bulk order
  return(furniture)
}




#' Make Metadata
#' 
#' Swipes the rat:group pairs from the big table in order to make the metadata
#'
#' @param data_tab 
#'
#' @return returns metadata as data frame
#' @export
#'
#' @examples
#' #Make metadata for either object (glut_tabs is properly populated)
#' metaData <- make_meta(glut_tabs)
make_meta <- function(data_tab){
  #This *should* work for another experiment as long as the big table has at 
  #least one cluster in it.
  ratG <- names(data_tab[2:9])
  
  #splits into list and makes storage for just the groups
  temp <- sapply(ratG, strsplit, ":")
  justG <- c()
  
  #Fills justG w/ the groups in the split
  for (i in 1:length(temp)){
    justG <- append(justG, temp[[i]][2]) 
  }
  
  #Makes data frame and labels it 
  met_fram <- data.frame(justG)
  rownames(met_fram) <- ratG
  colnames(met_fram) <- "group"  #<- Will likely need to be changed to a vector if sex also wanted
  
  #checks that the metadata is set up correctly
  cat("There should be a 'TRUE' under this line:\n")
  cat(as.character(all(rownames(met_fram) == colnames(counts)), "\n", sep = ""))
  cat("\n")
  
  #returns frame
  return(met_fram)
}


#' Get Pie Pieces
#' 
#' Returns chosen cluster's count table
#'
#' @param data_tab Aggregate table for the given object, returned by 
#'                 'get_cnt_tbls'
#' @param data_lst List of object's split cluster objects, returned by 
#'                 'make_cluster_list'
#' #@param clust_vect Vector of used clusters, returned by 'Kowalski'
#' @param object SE
#' @param dataname name of object, used a lot for once, if you put in the wrong
#'                 name, this function breaks
#'
#' @return the user's chosen count table
#' @export
#'
#' @examples
#' #Get count table for a cluster in the glut object
#' my_cnt_tbl <- prep_DEseq(glut_tab, glut, "Glut")
prep_DEseq <- function(data_tab, object, dataname){
  #prompt user for table
  cat("What cluster do you want to see a count table for?", 
      "Press 'v' to view options.", sep = "\n")
  piece = readline()
  
  #vectors to store cluster names
  data_lst <- c()
  clust_vect <- c()
  
  #populates data_lst 
  for (name in names(unique(data_tab))){
    if (str_detect(name, dataname)){
      data_lst <- append(data_lst, name)
    }
  }
  
  #splits data_lst into a list 
  temp <- sapply(data_lst, strsplit, dataname)
  
  #populates clust_vect using numbers from data_lst
  for (i in 1:(length(temp))){
    clust_vect <- append(clust_vect, temp[[i]][2]) 
  }  
  
  #removes possible NAs from clust_vect
  clust_vect <- clust_vect[!is.na(clust_vect)]  
  
  #Verifies 
  while (!((piece %in% data_lst) | (piece %in% clust_vect)) 
         & piece != "all"){
    if (piece != "v"){
      cat("Invalid input. Please enter a valid cluster or 'v'. \n") 
    }
    #Shows options
    cat("Your options are:\n")
    sapply(data_lst, print, quote = FALSE)
    sapply(clust_vect, print, quote = FALSE)
    print("all", quote = FALSE)
    print("v", quote = FALSE)
    #reruns prev line
    cat("What cluster do you want to see a count table for?", 
        "Press 'v' to view options.", sep = "\n")
    piece = readline()
  }
  
  #if a single num (like a 4), turns to something we can check for in the data Table
  if (piece %in% clust_vect){
    piece <- paste(dataname, piece, sep = "")
  }
  #same as prev but with all
  if (piece == "all"){
    piece <- paste("all", dataname, sep = "")
  }
  
  #Finds section
  tbl_sect <- match(piece, names(data_tab))
  
  #gets piece of the pie (table)
  tbl_piece <- data_tab[1:nrow(data_tab), (tbl_sect + 1):(tbl_sect + length(unique(object@meta.data[["ratID"]])))]
  return(tbl_piece)
}


#' Get Result Tibble
#' 
#' Creates result tables from count tables
#'
#' @param meta_tbl table of metadata
#' @param cnt_tbl count table you want a results table for 
#' @param object SE
#'
#' @return a result table
#' @export
#'
#' @examples
#' #Finds result tibble for count table glut0_cnts
#' my_res_tbl <- get_result_tbl(glut0_cnts, glut)
#' #Same as above but uses a different metadata table than the default
#' get_result_tbl(glut0_cnts, glut, my_meta_tbl)
get_result_tbl <- function(cnt_tbl, object, meta_tbl = metaData){
  #Gets user's input for contrast
  cat("Enter 2 groups, should be same as before")    #<- This may need to be changed if more wanted
  contrast_vect <- c("group", readline("Group1: "), readline("Group2: "))
  
  #verifies user input
  contrast_vect <- verify_factor(object, contrast_vect, "group")
  
  #Makes object
  cluster <- DESeqobjectFromMatrix(cnt_tbl,  #Data Table
                                    colData = meta_tbl,  #metaData
                                    design = ~ group)
  
  #Filters out significant parts
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
  
  #returns results
  return(results_tib)
}


#' Get Result Tibbles 2
#' 
#' This version combines prep_DEseq and get_result_tbl. Additionally, this has 
#' default parameters to make your life easier (half the time). The meta_tbl 
#' param can be left alone regardless of object unless you want to use different
#' metadata table. 
#'
#' @param data_tab Aggregate table for the given object, returned by 'get_cnt_tbls'
#' @param object SE
#' @param dataname name of object, used a lot for once
#' @param meta_tbl The table of metadata
#'
#' @return a result table
#' @export
#'
#' @examples
#' #makes result table for a cluster in the glut object
#' my_results <- run_DEseq()
#' #makes result table for a cluster in the GABA object
#' my_results <- run_DEseq(gaba_tabs, gaba, "Gaba")
#' #result table for a cluster in glut but w/ a different metadata table
#' my_results <- run_DEseq(meta_tbl = my_meta_tbl)
run_DEseq <- function(data_tab = glut_tabs, object = glut, 
                           dataname = "glut", meta_tbl = metaData){
  #Gets count table
  cnt_tbl <- prep_DEseq(data_tab, object, dataname)
  
  #gets and returns result table
  return(get_result_tbl(meta_tbl, cnt_tbl, object))
}


#' Make Mega Table
#'
#' @param object SE -- drakeeeee
#'
#' @return a 'MEGA' table of all of the chosen cluster's count tables
#' @export
#'
#' @examples
#' #makes an aggragate table combining all the count tables of user chosen 
#' #clusters
#' my_glut_cnts <- mk_mega_tbl(glut)
mk_mega_tbl <- function(object){
  
  #prompts user for what they want their data to be named
  cat("What would you like the name for your object to be?",
      "Good examples include:",
      "Glut",
      "GABA", sep = "\n")
  dataname = readline()
  
  #confirms the name
  print(paste("Are you sure you want your object to be named '", dataname, "'",
              "? (y/n)", sep = ""), quote = FALSE)
  confirm = readline()
  
  #gives user chance to change
  while (confirm != "y"){
    #invalid input
    if (confirm != "n"){
      confirm = readline("Invalid input. Please enter a 'y' or 'n'. ")
    }else {
      #prompts user again
      cat("What would you like the name for your object to be?",
          "Good examples include:",
          "Glut",
          "GABA", sep = "\n")
      dataname = readline()
      print(paste("Are you sure you want your object to be named '", dataname, "'",
                  "? (y/n)", sep = ""))
      confirm = readline()
    }
  }
  
  #Makes our data containers
  # clust_nams <- c()
  # clust_lst <- list()
  data_tabs <- matrix()
  
  #gets and stores the cluster variables that user defines
  #clust_nams <- Kowalski(object, dataname)
  
  #creates and stores seurot objects from split clusters
  #clust_lst <- make_clust_list(clust_nams, object, dataname)
  
  #makes the covetted MEGA table
  #data_tabs <- get_cnt_tbls(clust_lst, clust_nams, object, dataname)
  data_tabs <- get_cnt_tbls(make_clust_list(clust_nams, object, dataname),
                            Kowalski(object, dataname),
                            object,
                            dataname)
  
  #returns the MEGA table
  return(data_tabs)
}



#' Prep Pseudobulk
#' 
#' Takes in user input to create a list of Seurat objects with clusters chosen 
#' by the user.
#'
#' @param object Seurat object to perform analysis on 
#'
#' @return a list of Seurat objects with user chosen clusters
#' @export
#'
#' @examples
prep_psuedobulk <- function(object){
    #Prompts for threshold
    threshold = as.numeric(readline(
                  "What would you like the threshold for your object to be?: "))
    
    #validates input
    while (is.na(threshold)) {
        cat("Invalid input. Please enter a number.")
        threshold = as.numeric(readline(
                  "What would you like the threshold for your object to be?: "))
    }
    
    #runs filter_clusters and prints
    cat("Table:\n")
    
    #prompts user for factor
    factor = readline("What will your factor be?: ")
    
    #verifies that factor is somethning that exists within the object
    while (!(factor %in% colnames(object@meta.data))) {
        cat("Invalid input. Please try something like:\n",
            "group\n",    
            "sex\n", 
            "ratID\n", sep = "")
        factor = readline("What will your factor be?")
    }
    
    #stores the comparisons within the factor for ease of access
    factor_nams <- unique(object[[factor]])[, 1]
    cat("What would you like to compare?")
    cmprsn <- c(readline("Comparison1: "), 
                readline("Comparison2: "))
    
    #verifies user input
    cmprsn <- verify_factor(object, cmprsn, factor)
    
    #Makes base data table so we can import the data to the frame
    data_tbl <- table(object$seurat_clusters, object[[factor]][, 1])
    
    #Makes data frame to store
    clust_tab <- data.frame(matrix(ncol = length(factor_nams), 
                                   nrow = length(rownames(data_tbl))))
    
    #labels the frame's rows and columns
    colnames(clust_tab) <- factor_nams
    rownames(clust_tab) <- unique(levels(object@active.ident))
    
    #populates data frame
    for (i in 1:length(factor_nams)){
        clust_tab[, factor_nams[i]] <- data_tbl[, i]
    }
    
    #Very specific right now, Could probably be changed to reflect greater group 
    #size
    final <- clust_tab[clust_tab[, cmprsn[1]] > threshold &
                       clust_tab[, cmprsn[2]] > threshold, ] %>% 
                       select(cmprsn[1], cmprsn[2])
    print(final)
    cat("\n")
    
    #Makes vector of possible clusters to use
    poss <- rownames(final)
    
    #returns wanted clusters
    dVect <- c() #will be returned at function end
    
    #prompts user
    use_tbl = readline(paste("Would you like to keep the clusters in this",
                             "table? (y/n): "))
    #input validation
    while (use_tbl != 'y' & use_tbl != 'n') {
        use_tbl = readline("Invalid input. Please return a 'y' or 'n': ")
    }
    
    #yes, so sets dVect to poss
    if (use_tbl == 'y'){
      dVect <- poss
      
    #No, prompts user to ask for wanted clusters
    }else {
        #prompts user on cluster amt
        cs = as.numeric(readline("How many clusters would you like to use?: "))
        maximum <- length(poss)
        
        #
        while (is.na(cs)) {
            cat("Invalid input. Please enter a number.")
            cs = as.numeric(readline(paste("How many clusters would you like",
                                           "to use?: ")))
            cs <- as.numeric(cs)
        }
        
        #Basically means this is out of bounds
        while(maximum < cs | cs < 1){
            if (maximum < cs){
                cs = as.numeric(readline(paste("Choice exceeds current", 
                                               "clusters. Please use a",
                                               "lower number")))
            }else {
              cs = as.numeric(readline(paste("Choice is too small. Please use a",
                                             "greater number")))
            }
        }
        #max means all 
        if (cs == maximum){
            dVect <- poss
        }else {
            #prompts user for cs clusters by looping
            for (i in 1:cs){
                #Makes prompt
                prompt <- paste("Cluster ", i, ": ", sep = "")
                clust = readline(prompt)
                
                #loops until valid
                while (!(clust %in% poss) | (clust %in% dVect)){
                    cat("Invalid. Please choose a valid cluster that you have ",
                        "not already called.\n",
                        "Your options include:\n")
                    sapply(poss, print, quote = FALSE)
                    clust = readline(prompt)
                }
                
              #appends to dVect
              dVect <- append(dVect, clust)
            }
        }
    }
    
    #prints names so user can see
    print(dVect)
    
    #will be returned at function's end
    new_list <- list()
    
    #Makes cluster subsets derived from user input
    for (i in 1:length(dVect)) {
        
        #Makes val of new_list the subset
        new_list[[dVect[i]]] <- subset(object, idents = dVect[i])
        
        #Splitting stuff
        Idents(new_list[[dVect[i]]]) <- factor
        new_list[[dVect[i]]] <- subset(new_list[[dVect[i]]], idents = factor_nams)
    }
    
    #returns list
    return(new_list)
}

##---------------------------------------------------------------------------##


# Done for User # -----------------------------------------------------------
## Data Containers ---------------------------------------------------------


# Matrices to hold count tables
glut_tabs <- matrix()
gaba_tabs <- matrix()

# Will hold metadata, shouldn't need to change between objects
metaData <- data.frame()


## Tasks Done for User -------Makes MEGA Count Table ---------#---------
#glut data
# if (choice == 1 | choice == 2){
#   cat("Name your object something related to glut as this runs the code with glut\n")
#   glut_tabs <- mk_mega_tbl(glut)
# }
# #GABA data
# if (choice == 1 | choice == 3){
#   cat("Name your object something related to GABA as this runs the code with GABA\n")
#   gaba_tabs <- mk_mega_tbl(gaba)
# }
# 
# #metadata generated
# if (choice == 3){
#   metaData <- make_meta(gaba_tabs)
# }else {
#   metaData <- make_meta(glut_tabs)
# }

