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
                " two comparisons. \n", sep = "")
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
        #These lines will probably never trigger but it's good to have
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
