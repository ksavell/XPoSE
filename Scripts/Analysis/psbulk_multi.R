# The ps_bulk functions but modified to be able to take two factors

# By: Drake Thompson


# Info --------------------------------------------------------------------

#To be filled later


# Packages ----------------------------------------------------------------

# loads all required packages
library(dplyr)
library(purrr)
library(Seurat)
library(Libra)
library(here)
library(rlist)


# Data Loading ------------------------------------------------------------
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
#Loads data depending on user choice
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


# Functions ---------------------------------------------------------------

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
        while (!(vect[i] %in% valids)){
          #checks for case, converts to proper if not
          if (tolower(vect[i]) %in% valids_low){
            vect[i] <- valids[match(tolower(vect[i]), valids_low)]
          }else {
            #prompts user about 
            cat("Invalid input for Comparison", i, 
                ". Please check for typo.\n",
                "Valid inputs include:\n", sep = "")
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
      while (!(vect[i] %in% valids)){
        #checks for case 
        if (tolower(vect[i]) %in% valids_low){
          vect[i] <- valids[match(tolower(vect[i]), valids_low)]
        }else {
          cat("Invalid input for Comparison", i, 
              ". Please check for typo.\n",
              "Valid inputs include:\n", sep = "")
          sapply(valids, print, quote = FALSE)
          vect[i] = readline(paste("Comparsion", i, ": ", sep = ""))
        }
      }
    }
  }
  
  #Check to see if we need to rerun the function
  rerun <- FALSE
  
  #dupe check
  if (length(vect) > 1){
    for (i in (1 + mod):(length(vect) - 1)){
      while (vect[i] == vect[i + 1]){
        rerun <- TRUE
        cat("Dupe found at Comparsion", i + 1 - mod,
            ". Please use a different option within your set.\n",
            "Valid inputs include:\n", 
            "If you would like to delete this option, please enter 'd'. \n",
            sep = "")
        sapply(valids, print, quote = FALSE)
        cat("Please don't put in the same thing in...")
        vect[i + 1] = readline(paste("Comparsion", i + 1 - mod, ": ", sep = ""))
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


#' Verify Factors
#' 
#' General function that will take in a vector of user input and correct it using
#' a vector of the wanted factors' names. 
#'
#' @param object SE
#' @param vect the vector being verified, should have at least 2 elems
#' @param factors the factor vector to verify with. Within our object this will 
#'                contain chars from 'group' or 'sex'
#'
#' @return the user's (corrected) choice vector
#' @export
verify_factors <- function(object, vect, factors){
  #makes vector of valid options, normal and lowercase
  valids <- c()
  for (factor in factors){
    valids <- append(valids, unique(object[[factor]])[, 1])
  }
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
        while (!(vect[i] %in% valids)){
          #checks for case, converts to proper if not
          if (tolower(vect[i]) %in% valids_low){
            vect[i] <- valids[match(tolower(vect[i]), valids_low)]
          }else {
            #prompts user about 
            cat("Invalid input for Comparison", i, 
                ". Please check for typo.\n",
                "Valid inputs include:\n", sep = "")
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
      while (!(vect[i] %in% valids)){
        #checks for case 
        if (tolower(vect[i]) %in% valids_low){
          vect[i] <- valids[match(tolower(vect[i]), valids_low)]
        }else {
          cat("Invalid input for Comparison", i, 
              ". Please check for typo.\n",
              "Valid inputs include:\n", sep = "")
          sapply(valids, print, quote = FALSE)
          vect[i] = readline(paste("Comparsion", i, ": ", sep = ""))
        }
      }
    }
  }
  
  #dupe check
  if (length(vect) > 1){
    for (i in (1 + mod):(length(vect) - 1)){
      for (j in (2 + mod):length(vect)){
        while (vect[i] == vect[j]){
          cat("Dupe found at Comparsion", j - mod,
              ". Please use a different option within your set.\n",
              "Valid inputs include:\n", sep = "")
          sapply(valids, print, quote = FALSE)
          cat("Please don't put in the same thing in...")
          vect[j] = readline(paste("Comparsion", j - mod, ": ", sep = ""))
        }
      }
    }
  }
  
  #returns vect with changes (if any)
  return(vect)
}


#' Prep Pseudobulk
#' 
#' Takes in user input to create a list of Seurat objects with clusters chosen 
#' by the user. This version can take in up to 3 factors
#'
#' @param object Seurat object to perform analysis on 
#'
#' @return a list of Seurat objects with user chosen clusters
#' @export
#'
#' @examples
#' #make list of chosen clusters within glut object
#' clust_lst <- prep_psuedobulk(glut)
prep_multi_psbulk <- function(object){
  #Prompts for threshold
  cat("\n")
  threshold = as.numeric(readline(
    "What would you like the threshold for your object to be?: "))
  
  #validates input
  while (is.na(threshold)) {
    cat("Invalid input. Please enter a number.")
    threshold = as.numeric(readline(
      "What would you like the threshold for your object to be?: "))
  }
  
  #prompts for num of factors
  fact_num = as.numeric(readline(
    "How many factors would you like to use? (Max 3): "))
  
  #checks if user input is number
  while (is.na(fact_num)) {
    cat("Invalid input. Please enter a number.")
    fact_num = as.numeric(readline(
      "What would you like the threshold for your object to be?: "))
  }
  
  #ensures fact_num is valid size
  while (fact_num < 1 | fact_num > 3){
    cat("Invalid num. Please put in a number between 1 and 3. \n")
    
    #reprompts user
    fact_num = as.numeric(readline(
      "How many factors would you like to use? (Max 3): "))
    
    #checks if number again
    while (is.na(fact_num)) {
      cat("Invalid input. Please enter a number.")
      fact_num = as.numeric(readline(
        "What would you like the threshold for your object to be?: "))
    }
  }
  
  #will hold factor(s)
  factors <- c()
  
  #populates factor vect
  for (i in 1:fact_num){
  
    #prompts user for factor
    factor = readline(prompt = paste("What will your factor ", i,
                                     " be?: ", sep = ""))
    
    #verifies that factor is something that exists within the object
    while (tolower(factor) == tolower("ratID") |
           !(factor %in% colnames(object@meta.data)) |
           (tolower(factor) %in% tolower(factors))){
      #prevents use of ratID
      if (tolower(factor) == tolower("ratID")){
        cat("Biological replicates are required for analysis, so ratID cannot ",
            "be used as a factor. Please try something like:\n",
            "group\n",    
            "sex\n", 
            "What will factor ", i," be?: ", sep = "")
        factor = readline()
      #prevents dupes
      }else if (tolower(factor) %in% tolower(factors)){
        cat("Factor '", factor, "' already used. Please use another factor: ", sep = "")
        factor = readline()
      #makes this case insensitive
      }else if (tolower(factor) %in% tolower(colnames(object@meta.data))){
        factor <- colnames(object@meta.data)[
                  match(tolower(factor), tolower(colnames(object@meta.data)))]
      #User put in something invalid
      }else {      
        cat("Invalid input. Please try something like:\n",
          "group\n",    
          "sex\n", 
          "What will factor ", i," be?: ", sep = "")
        factor = readline()
      }
    }
    
    #Adds to vector
    factors <- append(factors, factor)
  }
  
  #will store things related to factors
  rat_list <- list()
  cmprsn <- list()
  
  #Index to track since we're in a for loop
  index <- 1
  
  #Gets comparisons for each factor
  for (factor in factors){
    #stores the comparisons within the factor for ease of access
    cmprsn[[factor]] <- c()
    #Makes sure length is good
    while (length(cmprsn[[factor]]) < 2){      
      cat("Please enter 2-3 comparisons for factor '", factor, "'.\n",
          "Press [enter] or 'q' on one of the prompts if you only want to make", 
          " two comparisons. \n",
          sep = "")
      cmprsn[[factor]] <- c(readline("Comparison1: "), 
                            readline("Comparison2: "),
                            readline("Comparison3 (optional): "))
      
      #removes if it's one of the chars below
      cmprsn[[factor]] <- cmprsn[[factor]][cmprsn[[factor]] != "" &
                                           cmprsn[[factor]] != "q" &
                                           cmprsn[[factor]] != " "]
    }
    
    #verifies user input
    cmprsn[[factor]] <- verify_factor(object, cmprsn[[factor]], factor)
    
    #"Makes" a table so that we can see how many rows it has. The num rows is 
    #stored, not the table
    rat_list <- list.append(rat_list,
                            factor = nrow(table(object$ratID,
                                                object[[factor]][, 1])[
                            table(object$ratID,
                                  object[[factor]][, 1])[,
                                          cmprsn[[factor]][1]] > 0 |
                            table(object$ratID,
                                  object[[factor]][, 1])[, 
                                          cmprsn[[factor]][2]] > 0, ]))
   
    #Confirmation it is the right name
    names(rat_list)[index] <- factor
    
    #increments index
    index <- index + 1
  }
  
  #Finds primary factor so data isn't wonky
  if (length(factors) > 1){
    #stores name of minimum
    minim <- c(names(rat_list)[1])
    
    #compares min to find primary factor
    for (i in 2:length(factors)){
      
      #New min found
      if (rat_list[[minim]] > rat_list[[i]]){
        minim <- c(names(rat_list)[i])
      
      #Tie between mins, user will be given choice
      } else if (rat_list[[minim]] == rat_list[[i]]){
        minim <- append(minim, names(rat_list)[i])
      }
    }
    
    #tie, user is given a choice
    if (length(minim) > 1){
      cat("No primary found. Please choose a factor that will help you filter ",
          "your clusters. \n",
          "Your options include:\n", sep = "")
      sapply(minim, print, quote = FALSE)
      factor = readline()
      
      while (!(factor %in% minim)) {
        cat("Invalid option. Please choose between the factors below:\n",
            sep = "")
        sapply(minim, print, quote = FALSE)
        factor = readline()
      }
    } else {
      factor <- minim[1]
    }
    
    #Changes elem so the right "hidden" factor(s) is stored
    names(rat_list)[match(factor, factors)] <- "min"
    
    #Reports to  user
    cat(factor, " chosen as primary factor. \n", sep = "")
  }else {
    #sets factor to only option
    factor <- factors[1]
  }
  
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
  if (length(cmprsn[[factor]]) == 2){
    final <- clust_tab[clust_tab[, cmprsn[[factor]][1]] > threshold &
                       clust_tab[, cmprsn[[factor]][2]] > threshold, ] %>% 
                       select(cmprsn[[factor]][1], cmprsn[[factor]][2])
  }else if (length(cmprsn[[factor]]) == 3){
    final <- clust_tab[(clust_tab[, cmprsn[[factor]][1]] > threshold &
                        clust_tab[, cmprsn[[factor]][2]] > threshold &
                        clust_tab[, cmprsn[[factor]][3]] > threshold), ] %>% 
                        select(cmprsn[[factor]][1], 
                               cmprsn[[factor]][2], 
                               cmprsn[[factor]][3])
  }else {
    cat("Something goofy has occurred. You have managed to make a table",
        "with either less than 2 or 4+ comparisons.\n",
        "Not sure how you managed that. Please try again.\n", sep = "")
    return()
  }
  print(final)
  cat("\n")
  
  #Makes vector of possible clusters to use
  poss <- rownames(final)
  
  #returns wanted clusters
  dVect <- c() 
  
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
    
  #No, so prompts user to ask for wanted clusters
  }else {
    #prompts user on cluster amt
    cs = as.numeric(readline("How many clusters would you like to use?: "))
    maximum <- length(poss)
    
    #Ensures a non numeric is not used
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
      
      #Ensures a non numeric is not used
      while (is.na(cs)) {
        cat("Invalid input. Please enter a number.")
        cs = as.numeric(readline(paste("How many clusters would you like",
                                       "to use?: ")))
        cs <- as.numeric(cs)
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
  
  #prints cluster names so user can see
  print(dVect)
  
  #will be returned at function's end, will contain the objs
  new_list <- list()
  
  #Makes cluster subsets derived from user input
  for (i in 1:length(dVect)) {
    
    #Makes val of new_list the subset
    new_list[[dVect[i]]] <- subset(object, idents = dVect[i])
    
    #Splitting stuff
    Idents(new_list[[dVect[i]]]) <- factor
    new_list[[dVect[i]]] <- subset(new_list[[dVect[i]]], idents = cmprsn[[factor]])
  }
  
  #Adds other factors as "hidden" data
  if (length(factors) > 1){
    new_list[["hidden"]] <- names(rat_list)[names(rat_list) != "min"]
  }
  
  #returns list
  return(new_list)
}



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
count_multi_psbulk <- function(data_lst, object){
    
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
    
    #transfers hidden 
    if (hidden){
      furniture[["hidden"]] <- data_lst[["hidden"]]
    }
    
    #returns list of count tables
    return(furniture)
}


#' Make Metadata
#' 
#' Swipes the rat:factor pairs from the big table in order to make the metadata
#'
#' @param tbl_lst list of the psuedobulk counts
#' @param object Seurat object to perform analysis on
#'
#' @return returns metadata as data frame
#' @export
#'
#' @examples
#' #Make metadata for the glut object
#' meta_data <- prep_DESeq(my_tbls, glut)
prep_mu_DESeq <- function(tbl_lst, object){
    
    #splits names into list of the ratID's and groups
    temp <- sapply(names(tbl_lst[[1]]), strsplit, ":")
    justF <- c()
    rat_nams <- c()
    
    #Fills justF w/ only the comparisons in the split
    for (i in 1:length(temp)){
        justF <- append(justF, temp[[i]][2]) 
        rat_nams <- append(rat_nams, temp[[i]][1])
    }
    
    #finds factor used in the list for user
    for (name in colnames(object@meta.data)){
        if (justF[1] %in% unique(object@meta.data[name][, 1])){
            factor <- name
        }
    }
    
    #Makes data frame and labels it 
    met_fram <- data.frame(justF)
    rownames(met_fram) <- names(tbl_lst[[1]])
    colnames(met_fram) <- factor  
    
    #Marks whether we have to transfer the vect in hidden or not
    hidden <- FALSE
    
    #makes new column for extra meta data
    if ("hidden" %in% names(tbl_lst)){
        #flags hidden as true
        hidden <- TRUE
      
        #loops through extra fectors
        for (i in 1:length(tbl_lst[["hidden"]])){
            #Makes vector to store what will be in new column
            new_col <- c()
            
            #loops through list of rat names get an index
            for (rat in rat_nams){
                #uses the index of the ratID to get what's within factor
                new_col <- append(new_col, 
                                  glut@meta.data[[tbl_lst[["hidden"]][i]]][
                                  match(rat, glut$ratID)][[1]])
            }
            
            #Adds vector as column and names it properly
            met_fram <- cbind(met_fram, factor = new_col)
            colnames(met_fram)[1 + i] <- tbl_lst[["hidden"]][i]
        }
      
        #Makes aggregate cols
        for (i in 1:length(tbl_lst[["hidden"]])){
            #Will hold aggregates
            new_col <- c()
            
            #makes the data within
            for (j in 1:nrow(met_fram)){
                #makes factor-hidd pairs
                new_col <- append(new_col, 
                                  paste(met_fram[, factor][j], "_", 
                                        met_fram[, tbl_lst[["hidden"]][i]][j], 
                                        sep = ""))
            }
            #Adds vector as column and names it properly
            met_fram <- cbind(met_fram, factor = new_col)
            colnames(met_fram)[1 + 
                               length(tbl_lst[["hidden"]]) + i] <- paste(factor,
                                          "_", tbl_lst[["hidden"]][i], sep = "")
        }
      
        #for tri-factor sets
        if (length(tbl_lst[["hidden"]]) == 2){
            
            ind <- 1
            #targets aggregate vectors
            for (agg in colnames(met_fram)[4:5]){
                
                #Will hold tri-factor aggregates
                new_col <- c()
                
                #swipes the datapiece not in the aggregate 
                singl <- colnames(met_fram)[1:3][
                          !(colnames(met_fram)[1:3] %in% str_split_1(agg, "_"))]
                
                #makes the data within
                for (j in 1:nrow(met_fram)){
                    #makes tri-factor data
                    new_col <- append(new_col, 
                                      paste(met_fram[, agg][j], "_", 
                                            met_fram[, singl][j], 
                                            sep = ""))
                }
                
                #Adds vector as column and names it properly
                met_fram <- cbind(met_fram, factor = new_col)
                colnames(met_fram)[5 + ind] <- paste(agg, "_", singl, sep = "")
                
                #increments ind
                ind <- ind + 1
            }
        }
    }
    
    #checks that the metadata is set up correctly
    cat("There should be at least one 'TRUE' on each line:\n")
    cat(as.character(unique(met_fram[,1]) %in% 
                     unique(object@meta.data[factor][, 1])), "\n")
    #Checks possible hidden columns
    if (hidden){
        for (i in 1:length(tbl_lst[["hidden"]])) {
            cat(as.character(unique(met_fram[, 1 + i]) %in% 
                             unique(object@meta.data[
                             tbl_lst[["hidden"]][i]][, 1])), "\n")
        }
    }
    cat("\n")
    
    #returns frame
    return(met_fram)
}


#You are here. Helllo tomorrow me, hope you slept goood last night


#' Run DEseq
#'
#' @param tbl_lst list of the psuedobulk counts
#' @param meta_tbl table of metadata
#' @param object Seurat object to perform analysis on
#'
#' @return returns aggegate list of the result tables (they are tibbles)
#' @export
#'
#' @examples
#' #Makes result tibbles for my_tbls list
#' my_tibbs <- run_DEseq(my_tbls, meta_data, glut)
run_DEseq <- function(tbl_lst, meta_tbl, object){
  #our base list to hold tibbles, will be returned at function's end
  agg_res <- list()
  
  #finds factor used in the list for user
  for (name in colnames(object@meta.data)){
    if (unique(meta_tbl[,1])[1] %in% unique(object@meta.data[name][, 1])){
      factor <- name
    }
  }
  
  #Gets user's input for contrast
  cat("Enter 2 comparisons, should be same as before")  #<- This may need to 
  #   be changed if more 
  #   wanted
  contrast_vect <- c(factor, readline("Experimental Group: "), 
                     readline("Control Group: "))
  
  
  #Makes factor a formula for the design
  form_fact <- paste("~", factor, sep = "")
  form_fact <- as.formula(form_fact)
  
  #verifies user input
  contrast_vect <- verify_factor(object, contrast_vect, factor)
  
  #Makes result table for each count table in list
  for (i in 1:length(tbl_lst)){
    #Makes object
    cluster <- DESeqDataSetFromMatrix(tbl_lst[[i]],  #Data Table
                                      colData = meta_tbl,  #metaData
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
    
    #returns results
    agg_res <- list.append(agg_res, clust_nam = results_tib)
    names(agg_res)[i] <- names(tbl_lst)[i]
  }
  
  return(agg_res)
}
