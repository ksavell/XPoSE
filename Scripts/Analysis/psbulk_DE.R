# Remakes the remade base file with code sections, minor improvements and such

# By: Drake Thompson, Padmashri Saravanan, Katherine Savell

# Info --------------------------------------------------------------------

# This code will bringyou through the main steps of the initial differential
# analysis. This can be used to make plots or just look at the data differently.
# There are 4 functions with another that just helps some of them. 
# The functions should be done in this order:
# 1. prep_pseudobulk
# 2. count_psuedobulk
# 3. prep_DEseq
# 4. run_DEseq
# The function after usually needs a material acquired from the previous. 


# Loading # -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
# loads all required packages
library(dplyr)
library(purrr)
library(Seurat)
library(Libra)
library(here)
library(rlist)


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


#Functions -----------------------------------


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
#' by the user.
#'
#' @param object Seurat object to perform analysis on 
#'
#' @return a list of Seurat objects with user chosen clusters
#' @export
#'
#' @examples
#' #make list of chosen clusters within glut object
#' clust_lst <- prep_psuedobulk(glut)
prep_psuedobulk <- function(object){
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
    
    #prompts user for factor
    factor = readline("What will your factor be?: ")
    
    #verifies that factor is somethning that exists within the object
    while (!(factor %in% colnames(object@meta.data))) {
        cat("Invalid input. Please try something like:\n",
            "group\n",    
            "sex\n", 
            "ratID\n", sep = "")
        factor = readline("What will your factor be? ")
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
    
    #Prints table for user pleasure
    cat("\nTable:\n")
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
      
    #No, so prompts user to ask for wanted clusters
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
        new_list[[dVect[i]]] <- subset(new_list[[dVect[i]]], idents = cmprsn)
    }
    
    #returns list
    return(new_list)
}


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
prep_DESeq <- function(tbl_lst, object){
    
    #splits names into list of the ratID's and groups
    temp <- sapply(names(tbl_lst[[1]]), strsplit, ":")
    justF <- c()
    
    #Fills justF w/ only the comparisons in the split
    for (i in 1:length(temp)){
        justF <- append(justF, temp[[i]][2]) 
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
    colnames(met_fram) <- factor  #<- Will likely need to be changed if multiple
                                  #   factors wanted
    
    #checks that the metadata is set up correctly
    cat("There should be at least one 'TRUE' under this line:\n")
    cat(as.character(unique(met_fram[,1]) %in% unique(object@meta.data[factor][, 1])), "\n")
    cat("\n")
    
    #returns frame
    return(met_fram)
}


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
    
    #list of possible "formulas", will need to be updated eventually
    poss_des <- list(group = ~group, 
                     sex = ~sex, 
                     ratID = ~ratID, 
                     population = ~population,
                     experience = ~experience, 
                     region = ~region)
   
    #match factor w/ formula
    for (f in names(poss_des)){
        if (f == factor){
            des <- poss_des[[f]]
        }
    }
    
    #verifies user input
    contrast_vect <- verify_factor(object, contrast_vect, factor)
    
    #Makes result table for each count table in list
    for (i in 1:length(tbl_lst)){
        #Makes object
        cluster <- DESeqDataSetFromMatrix(tbl_lst[[i]],  #Data Table
                                          colData = meta_tbl,  #metaData
                                          design = des)

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


