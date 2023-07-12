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
