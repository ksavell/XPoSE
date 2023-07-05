#Remakes the remade base file with code sections, minor improvements and such
#By: Drake Thompson, Padma Saravanan, Katherine Savell

# Info --------------------------------------------------------------------

#This doesn't make an Upset plot but sets up a lot of the steps. Everything up 
#to the last 3 functions is done for you. 
#This makes two main things: 
#           a 'MEGA' table that can divided into subtables 
#                       (use get_pie_piece to get subtables)
#           metadata table (stored  as metaData) for making result tables
#                       (use get_result_tbl or get_result_2bl)
#Making individual count and result tables will have to be done by user
#This is done w/ the three remaining functions:
#                                              get_pie_piece
#                                              get_result_tbl
#                                              get_result_2bl

#Note: This code is broken into various sections. You can collapse sections and 
#look at specific chunks of code. There is also a Jump To Command that takes you 
#to the specific sections under 'Code' on the top menu bar.


# Loading # -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
#loads all required packages
library(dplyr)
library(purrr)
library(magrittr)
library(tester)
library(Matrix)
library(lmtest)
library(tidyselect)
library(DESeq2)
library(Seurat)
library(blme)
library(glmmTMB)
library(lme4)
library(lmerTest)
library(matrixStats)
library(methods)
library(stats)
library(Rdpack)
library(Libra)
library(pheatmap)
library(tibble)
library(UpSetR)
library(ComplexHeatmap)
library(ComplexUpset)
library(data.table)


## Data Loading ------------------------------------------------------------
#Loads data depending on user choice
prompt <- paste("Use both data sets?\n", 
                "1. Yes\n",
                "2. Just Glut\n",
                "3. Just GABA\n")
cat(prompt)
choice = readline()

#validates input, loops if invalid
while (!(choice %in% 1:3)) {
  #"error" message
  print("Invalid Input. Please return a 1, 2, or 3")
  
  #asks fo input again
  cat(prompt)
  choice = readline()
}

#Loads data, you may need to change below code depending on device. 
##-----------------------------------------------##
#Loads data depending on user choice
if (choice == 1 | choice == 2){
  to_load = readline("Load glut? (y/n): ")
  if (to_load == "y"){
    load("data/glut_subset_forthpipe_11102022.RData")
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
#SE refers to 'self explanatory'
#Functions are collapsable if you want to just read descriptions and this section
#feels cluttered.

## List of Functions that May Require Changes ------------------------------
##------------------------##
#List of functions that may need to be changed with a new experiment
##------------------------##
###verify_factor:
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
###get_pie_piece: 
#              The basis of this function findind the piece is on ratID so if its
#              location changes in the dataset, problems may occur. Please be aware
#              of that!





##Functions -----------------------------------


#' Verify Factor
#' 
#' General function that will take in a vector of user input and correct it using
#' a vector of the wanted factor's names. 
#'
#' @param dataset SE
#' @param vect the vector being verified
#' @param factor the factor to verify with. Within our dataset this will be 
#'               'group' or 'sex'
#'
#' @return the user's (corrected) choice vector
#' @export
verify_factor <- function(dataset, vect, factor){
  
  #makes vector of valid options
  valids <- unique(dataset[[factor]])[,1]
  
  #Loops vector
  for (i in 1:length(vect)){
    #There is a case where group is the first elem so this ensure output to user
    #is consistent
    mod <- 0 
    if (i == 1){
      if (vect[i] != factor){
        #Loops until valid input
        while (!(vect[i] %in% valids)){
          cat("Invalid input for Group", i, ". Please check for typo.\n",
              "Valid inputs include:\n", sep = "")
          sapply(valids, print, quote = FALSE)
          vect[i] = readline(paste("Group", i, ": ", sep = ""))
        }
      #If group = 1, i needs subtraction
      } else {
        mod <- 1
      }
    }else {
      #Loops until valid input
      while (!(vect[i] %in% valids)){
        cat("Invalid input for Group", (i + mod), ". Please check for typo.\n",
            "Valid inputs include:\n", sep = "")
        sapply(valids, print, quote = FALSE)
        vect[i] = readline(paste("Group", (i + mod), ": ", sep = ""))
      }
    }
  }
  
  #returns vect with changes (if any)
  return(vect)
}


#' Filter Clusters
#' 
#' Makes table of data filtered by wanted THRESHOLD and comparison
#'
#' @param dataset SE
#' @param threshold min amt required to be used. Decided by user
#'
#' @return the filtered table
#' @export
#'
#' @examples
#' #Makes table of glut clusters w/ counts > 100
#' filter_clusters(glut, 100)
#' #Makes table of GABA clusters w/ counts > 50
#' filter_clusters(gaba, 50)
filter_clusters <- function(dataset, threshold){
  #variables to reduce complexity
  grp_nams <- unique(dataset$group)
  cmprsn <- c(readline("Group1: "), #I could require inputs to be from grp_nams
              readline("Group2: "))
  
  #verifies user input
  cmprsn <- verify_factor(dataset, cmprsn, "group")
  
  #Makes base data table so we can import the data to the frame
  data_tbl <- table(dataset$seurat_clusters, dataset$group)#May change between datasets
  
  #Makes data frame to store
  df <- data.frame(matrix(ncol=length(grp_nams), nrow=length(rownames(data_tbl))))
  
  #labels the frame's rows and columns
  colnames(df) <- grp_nams
  rownames(df) <- rownames(data_tbl)
  
  #populates data frame
  for (i in 1:length(grp_nams)){
    df[,grp_nams[i]] <- data_tbl[,i]
  }
  
  #Very specific right now, Could probably be changes to reflect greater group size
  final <- df[df[,cmprsn[1]] > threshold & df[,cmprsn[2]] > threshold,] %>% 
    select(cmprsn[1], cmprsn[2])
  
  #Returns wanted comparison
  return(final)
}


#' Analysis
#' 
#' Helper function for Kowalski that gets the wanted threshold and resulting 
#' table from filter_clusters. 
#'
#' @param dataset SE
#' @param dataname name of the dataset being used
#'
#' @return table made from filter_clusters
#' @export
analysis <- function(dataset, dataname){
  
  #Prompts for threshold
  prompt <- paste("What would you like the threshold for your", dataname,
                  "data to be?: " )
  thres = readline(prompt)
  thres <- as.numeric(thres) #Code will break without this line
  
  #runs filter_clusters and prints
  cat("Table:\n")
  temp <- filter_clusters(dataset, thres)
  print(temp)
  cat("\n")
  
  #returns temp so it can be stored
  return(temp)
}



#' Fill Vect
#' 
#' Fills the respective vector with the clusters of user's choosing. 
#'
#' @param possible The possible clusters to choose from, derived from data frame
#'                 created by filter_clusters
#'
#' @return the vector of chosen clusters
#' @export
fill_vect <- function(possible){
  dVect <- c() #will be returned at function end
  
  #calls helper function to prompt
  use_tbl = readline("Would you like to keep the clusters in this table? (y/n): ")
  #input validation
  while (use_tbl != 'y' & use_tbl != 'n') {
    use_tbl = readline("Invalid input. Please return a 'y' or 'n': ")
  }
  
  #yes, so sets dVect to poss
  if (use_tbl == 'y'){
    dVect <- possible
    
  #No, prompts user to ask for wanted clusters
  }else {
    #prompts user on cluster amt
    cs = readline("How many clusters would you like to use?: ")
    max <- length(possible)
    
    #Basically means this is out of bounds
    while(max < cs | cs < 1){
      cs = readline("Choice exceeds current clusters. Please use a lower number")
    }
    #Calls self to restart, max means all 
    if (cs == max){
      #reruns function
      dVect <- fill_vect(possible)
      
    }else {
      #prompts user for cs clusters by looping
      for (i in 1:cs){
        #Makes prompt
        prompt <- paste("Cluster ", i, ": ", sep = "")
        clust = readline(prompt)
        
        #loops until valid
        while (!(clust %in% possible) | (clust %in% dVect)){
          cat("Invalid. Please choose a valid cluster that you have not already called.")
          clust = readline(prompt)
        }
        
        #appends to dVect
        dVect <- append(dVect, clust)
      }
    }
  }
  
  #returns resulting vector
  return(dVect)
}


#' Kowalski
#' "Decides" which clusters to use for each dataset. 
#' Uses analysis and fill_vect
#' 
#'
#' @param dataset SE
#' @param dataname The name of the dataset being used. Funneled to 'analysis'
#'
#' @return a vector of the cluster "names"
#' @export
Kowalski <- function(dataset, dataname){
  
  #Calls other functions to generate table
  #data_tut <- analysis(dataset, dataname)
  #Prompts for threshold
  prompt <- paste("What would you like the threshold for your", dataname,
                  "data to be?: " )
  thres = readline(prompt)
  thres <- as.numeric(thres) #Code will break without this line
  
  #runs filter_clusters and prints
  cat("Table:\n")
  temp <- filter_clusters(dataset, thres)
  print(temp)
  cat("\n")
  
  # #returns temp so it can be stored
  # return(temp)
  
  #Makes vector of possible clusters to use
  #poss <- rownames(data_tut)
  poss <- rownames(temp)
  
  #returns wanted clusters
  return(fill_vect(poss))
}


#' Make Cluster List
#'
#' @param clust_vect vector containing the clusters chosen in Kowalski
#' @param dataset SE
#' @param dataname The name of the dataset being used ("Glut" or "GABA")
#'
#' @return list of clusters (Seurot objects)
#' @export
make_clust_list <- function(clust_vect, dataset, dataname){
  #will be returned at function's end
  new_list <- list()
  #What the groups will be, user choice
  cat("Enter 2 groups for idents, should be same as before")    #<- This may need to be changed if more wanted
  filtration <- c(readline("Group1: "), readline("Group2: "))
  
  filtration <- verify_factor(dataset, filtration, "group")
  
  #Makes cluster subsets derived from user input
  for (i in 1:length(clust_vect)) {
    #converts cluster to char just in case
    clust_vect[i] <- as.character(clust_vect[i])
    
    #Makes temp for the name as it's frequently used
    tempNme <- paste(dataname, clust_vect[i], sep = "")
    
    #Makes val of new_list the subset
    new_list[[tempNme]] <- subset(dataset, idents = clust_vect[i])
    
    #Splitting stuff
    Idents(new_list[[tempNme]]) <- "group"
    new_list[[tempNme]] <- subset(new_list[[tempNme]], idents = filtration)
  }
  
  #returns list
  return(new_list)
}



#' Get Count Tables
#'
#' @param data_lst list of searot objects for dataser, generated by 
#'                 'make_cluster_list'
#' @param clust_vect vector of clusters for dataset, generated by 'Kowalski'
#' @param dataset SE, used to make allset
#' @param dataname name of data, used to make all_tbl
#'
#' @return the "MEGA" table that contains all the count tables
#' @export
#'
#' @examples
#' #Makes count tables for the requested clusters for glut dataset
#' mega_tbl <- get_cnt_tbls(glut_lst, glut_dat, glut, "Glut")
get_cnt_tbls <- function(data_lst, clust_vect, dataset, dataname){
  
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
  allset <- FindClusters(dataset, resolution = 0)
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
#If this function is too long for you, I can separate it into 2


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
#' #Make metadata for either dataset (glut_tabs is properly populated)
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
#' @param data_tab Aggregate table for the given dataset, returned by 
#'                 'get_cnt_tbls'
#' @param data_lst List of dataset's split cluster objects, returned by 
#'                 'make_cluster_list'
#' #@param clust_vect Vector of used clusters, returned by 'Kowalski'
#' @param dataset SE
#' @param dataname name of dataset, used a lot for once, if you put in the wrong
#'                 name, this function breaks
#'
#' @return the user's chosen count table
#' @export
#'
#' @examples
#' #Get count table for a cluster in the glut dataset
#' my_cnt_tbl <- get_pie_piece(glut_tab, glut, "Glut")
get_pie_piece <- function(data_tab, dataset, dataname){
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
  tbl_piece <- data_tab[1:nrow(data_tab), (tbl_sect + 1):(tbl_sect + length(unique(dataset@meta.data[["ratID"]])))]
  return(tbl_piece)
}


#' Get Result Tibble
#' 
#' Creates result tables from count tables
#'
#' @param meta_tbl table of metadata
#' @param cnt_tbl count table you want a results table for 
#' @param dataset SE
#'
#' @return a result table
#' @export
#'
#' @examples
#' #Finds result tibble for count table glut0_cnts
#' my_res_tbl <- get_result_tbl(glut0_cnts, glut)
#' #Same as above but uses a different metadata table than the default
#' get_result_tbl(glut0_cnts, glut, my_meta_tbl)
get_result_tbl <- function(cnt_tbl, dataset, meta_tbl = metaData){
  #Gets user's input for contrast
  cat("Enter 2 groups, should be same as before")    #<- This may need to be changed if more wanted
  contrast_vect <- c("group", readline("Group1: "), readline("Group2: "))
  
  #verifies user input
  contrast_vect <- verify_factor(dataset, contrast_vect, "group")
  
  #Makes object
  cluster <- DESeqDataSetFromMatrix(cnt_tbl,  #Data Table
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
#' This version combines get_pie_piece and get_result_tbl. Additionally, this has 
#' default parameters to make your life easier (half the time). The meta_tbl 
#' param can be left alone regardless of dataset unless you want to use different
#' metadata table. 
#'
#' @param data_tab Aggregate table for the given dataset, returned by 'get_cnt_tbls'
#' @param dataset SE
#' @param dataname name of dataset, used a lot for once
#' @param meta_tbl The table of metadata
#'
#' @return a result table
#' @export
#'
#' @examples
#' #makes result table for a cluster in the glut dataset
#' my_results <- get_result_2bl()
#' #makes result table for a cluster in the GABA dataset
#' my_results <- get_result_2bl(gaba_tabs, gaba, "Gaba")
#' #result table for a cluster in glut but w/ a different metadata table
#' my_results <- get_result_2bl(meta_tbl = my_meta_tbl)
get_result_2bl <- function(data_tab = glut_tabs, dataset = glut, 
                           dataname = "glut", meta_tbl = metaData){
  #Gets count table
  cnt_tbl <- get_pie_piece(data_tab, dataset, dataname)
  
  #gets and returns result table
  return(get_result_tbl(meta_tbl, cnt_tbl, dataset))
}


#' Make Mega Table
#'
#' @param dataset SE
#'
#' @return a 'MEGA' table of all of the chosen cluster's count tables
#' @export
#'
#' @examples
#' #makes an aggragate table combining all the count tables of user chosen 
#' #clusters
#' my_glut_cnts <- mk_mega_tbl(glut)
mk_mega_tbl <- function(dataset){
  #prompts user for what they want their data to be named
  cat("What would you like the name for your dataset to be?",
      "Good examples include:",
      "Glut",
      "GABA", sep = "\n")
  dataname = readline()
  
  #confirms the name
  print(paste("Are you sure you want your dataset to be named '", dataname, "'",
              "? (y/n)", sep = ""), quote = FALSE)
  confirm = readline()
  
  #gives user chance to change
  while (confirm != "y"){
    #invalid input
    if (confirm != "n"){
      confirm = readline("Invalid input. Please enter a 'y' or 'n'. ")
    }else {
      #prompts user again
      cat("What would you like the name for your dataset to be?",
          "Good examples include:",
          "Glut",
          "GABA", sep = "\n")
      dataname = readline()
      print(paste("Are you sure you want your dataset to be named '", dataname, "'",
                  "? (y/n)", sep = ""))
      confirm = readline()
    }
  }
  
  #Makes our data containers
  clust_nams <- c()
  clust_lst <- list()
  data_tabs <- matrix()
  
  #gets and stores the cluster variables that user defines
  clust_nams <- Kowalski(dataset, dataname)
  
  #creates and stores seurot objects from split clusters
  clust_lst <- make_clust_list(clust_nams, dataset, dataname)
  
  #makes the covetted MEGA table
  data_tabs <- get_cnt_tbls(clust_lst, clust_nams, dataset, dataname)
  
  #returns the MEGA table
  return(data_tabs)
}



##---------------------------------------------------------------------------##


# Done for User # -----------------------------------------------------------
## Data Containers ---------------------------------------------------------


#Matrices to hold count tables
glut_tabs <- matrix()
gaba_tabs <- matrix()

#Will hold metadata, shouldn't need to change between datasets
metaData <- data.frame()


## Tasks Done for User -------Makes MEGA Count Table ---------#---------
#glut data
if (choice == 1 | choice == 2){
  cat("Name your dataset something related to glut as this runs the code with glut\n")
  glut_tabs <- mk_mega_tbl(glut)
}
#GABA data
if (choice == 1 | choice == 3){
  cat("Name your dataset something related to GABA as this runs the code with GABA\n")
  gaba_tabs <- mk_mega_tbl(gaba)
}

#metadata generated
if (choice == 3){
  metaData <- make_meta(gaba_tabs)
}else {
  metaData <- make_meta(glut_tabs)
}


