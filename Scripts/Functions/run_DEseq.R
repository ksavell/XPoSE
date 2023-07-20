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
    source("~/Documents/GitHub/XPoSE/Scripts/Functions/verify_factor.R")
    
    #our base list to hold tibbles, will be returned at function's end
    agg_res <- list()
    
    #factor to use in contrast
    fctor <- ""
    
    #Users can append to hidden. This code prevents more than 2 factors
    #from existing in hidden when it is used. Also if they try to use ratID. 
    if ("hidden" %in% names(tbl_lst)){
        #prevents ratID use
        tbl_lst[["hidden"]] <- tbl_lst[["hidden"]][tbl_lst[["hidden"]] != "ratID"]
        
        #ensure only 2 extra factors
        if (length(tbl_lst[["hidden"]]) > 2){
            tbl_lst[["hidden"]] <- tbl_lst[["hidden"]][1:2]
        }  
    }
    
    #prompt user for used factor
    while (!(fctor %in% colnames(meta_tbl))){
        #invalid input
        if (fctor != ""){
            cat("Invalid input. Please choose one of the options below:\n",
                sep = "")
        #first hit/they entered nothing
        }else {
             cat("Please choose a factor or factors to do DESeq on. ",
                 "Options are below:\n", sep = "")
        }
        #displays options and allows input
        sapply(colnames(meta_tbl), print, quote = FALSE)
        cat("Underscores denote a combination of factors, please include ",
        "them if you're trying to use multiple factors.\n", sep = "")
        fctor = readline()
    }
    
    #Makes factor a formula for the design
    form_fact <- paste("~", fctor, sep = "")
    form_fact <- as.formula(form_fact)
    
    #Gets user's input for contrast
    cat("Enter 2 comparisons from your factor. ",
        "Comparisons should be entered in same order as before. ",
        "Examples below.\n")
    sapply(unique(meta_tbl[, fctor]), print, quote = FALSE)
    contrast_vect <- c(fctor, readline("Experimental Group: "),
                              readline("Control Group: "))
    
    #verifies user input
    if (fctor %in% colnames(meta_tbl)[1:(length(tbl_lst[["hidden"]]) + 1)]){
        #factor exists in metadata and can verified traditionally
        contrast_vect <- verify_factor(object, contrast_vect, fctor)
        
    #uses composite, gotta check with the list
    } else {
        #runs through every part except the factor
        for (i in 2:3){
            #validates input
            while (!(contrast_vect[i] %in% unique(meta_tbl[, fctor])) | 
                    (i == 3 & contrast_vect[2] == contrast_vect[3])) {
                #dupe
                if (i == 3 & contrast_vect[2] == contrast_vect[3]){
                    cat("Dupe found at Comparison 3. Please input something", 
                        "different from the list below:\n", sep = "")
                #invalid
                }else {
                    cat("Invalid comparison used. Please input something from the ",
                        "list below:\n", sep = "")
                }
                
                #shows options and gives tip before reprompt
                sapply(unique(meta_tbl[, fctor]), print, quote = FALSE)
                cat("Don't forget to include underscores!\n")
                
                #reprompts user
                contrast_vect[i] = readline(paste("Comparison", i - 1, ": ", sep = ""))
            }
        }
    }
    
    #Makes result table for each count table in list (excludes hidden)
    for (i in 1:length(tbl_lst)){
        #ensures hidden isn't touched
        if (names(tbl_lst)[i] != "hidden"){
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
            
            #appends to list of tibbles
            agg_res <- list.append(agg_res, clust_nam = results_tib)
            names(agg_res)[i] <- names(tbl_lst)[i]
        }
    }
    #returns results
    return(agg_res)
}
