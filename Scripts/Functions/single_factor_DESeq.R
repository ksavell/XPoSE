#' Single Factor DeSeq
#' 
#' This takes a single cluster and a single factor within a Seurat object and 
#' creates both a DeSeq and results table for the object for further analysis
#'
#' @param object The Seurot object to do analysis on
#' @param comp_vect a vector holding your factor and levels within the factor.
#'                  The first element of the vector should be your factor and 
#'                  the subsequent elements should be levels within you want to
#'                  compare.
#'                  
#'                  For example, if you wanted to analyze differences between 
#'                  Active and Non-active within your object, this vector would
#'                  look like:
#'                  
#'                  comp_vect <- c("group", "Active", "Non-active")
#'                  
#' @param cluster the cluster within your object you want to analyze
#' @param min_cell mimimum count that a cluster of each rat must have for that 
#'                 cluster to be considered in analysis. Has a minimum value of 
#'                 0.
#' @param min_rat A minimum number of rats you would want for your comparison.
#'                Cannot be lower than 2. Has a minimum value of 2.
#'
#' @return a list containing the DESeq object and the results table
#' @export
#'
#' @examples
single_factor_DESeq <- function(object, comp_vect, cluster, min_cell = 10, 
                                min_rat = 2, keep_dds = FALSE){
    library(Seurat)
    library(Libra)
    library(dplyr) 
    library(tibble)
    library(DESeq2)
    library(rlang)
    library(rlist)
    
    #checks that cluster exists
    if (!(cluster %in% unique(object$cluster_name))){
        stop("Cluster '", deparse(substitute(cluster)), "' not found within ",
             "object: ", deparse(substitute(object)), ".")
    }
    
    #ensures min_rat is at least 2
    if (min_rat < 2){
        min_rat <- 2
    }
    #esnures min_cell is at least 0
    if (min_cell < 0){
        min_cell <- 0
    }
    
    if (!(comp_vect[1] %in% colnames(object@meta.data))){
        stop("Factor, ", comp_vect[1], ", not found within object '", 
             deparse(substitute(object)), "'.")
    }
    if (FALSE %in% (comp_vect[2:3] %in% object[[comp_vect[1]]][, 1])){
        stop("Vector '", deparse(substitute(comp_vect)), 
             "' contains comparisons not in factor '", comp_vect[1], "'.")
    }
    
    #subsets objects by the comparisons
    Idents(object) <- "cluster_name"
    sub_obj <- subset(object, idents = cluster)
    Idents(sub_obj) <- comp_vect[1]
    sub_obj <- subset(sub_obj, idents = comp_vect[2:3])
    
    #makes table to display rats that will be included as well as their counts
    #and characteristics
    t_tbl <- data.frame(matrix(nrow = 1, ncol = 4))
    
    #table to pull relevant rats from 
    rat_g_c <- table(sub_obj[[comp_vect[1]]][, 1],  
                     sub_obj$cluster_name, sub_obj$ratID) 
    #names 
    colnames(t_tbl) <- c(comp_vect[1], "ratID", "counts", 
                         "Included")
    
    #the row of the table
    curr_row <- 1
    #loops through our rats and comparisons to fill table
    for (rat in unique(sub_obj$ratID)){
        for (level in comp_vect[2:3]){
            if (rat_g_c[level, ,rat] != 0){
                #data we fill row with
                t_tbl[curr_row, ] <- list(level, rat, 
                                          rat_g_c[level, ,rat], 
                                          rat_g_c[level, ,rat] > 
                                              min_cell)
                
                #updates name and index
                rownames(t_tbl)[curr_row] <- paste(rat, level, sep = ":")
                
                #updates curr_row
                curr_row <- curr_row + 1
            }
        }
    }
    
    #makes table to display holding the number of rats within our sample that 
    #have characteristics of a comparison over the total number of rats in the 
    #data provided with those same characteristics
    incl_tbl <- data.frame(matrix(0, nrow = 1, ncol = 2))
    rownames(incl_tbl) <- comp_vect[1]
    colnames(incl_tbl) <- comp_vect[2:3]
    
    rat_ttls <- c()
    #draws data from the original table
    for (i in 1:nrow(t_tbl)){
        #draws data from original table
        if (t_tbl[i, "Included"]){
            incl_tbl[, t_tbl[i, comp_vect[1]]] <- 
                                    incl_tbl[, t_tbl[i, comp_vect[1]]] + 1
        }
        
        #finds totals
        if (!(t_tbl[i, comp_vect[1]] %in% names(rat_ttls))){
            rat_ttls[t_tbl[i, comp_vect[1]]] <- 1
        }else {
            rat_ttls[t_tbl[i, comp_vect[1]]] <- 
                rat_ttls[t_tbl[i, comp_vect[1]]] + 1
        }
        
    }
    
    #flag is here as we want to print the warning after printing the table
    #but counts are easier to check prior to the step that adds the / total
    low_samp_flag <- FALSE
    #loops through table to find with min_rat
    for (level in comp_vect[2:3]){
            if (incl_tbl[,level] < min_rat){
                low_samp_flag <- TRUE
            }
            if (incl_tbl[,level] == 0){
               stop("Level, ", level, ", within factor ", comp_vect[1], 
                    " was fully excluded with chosen min-cell. Please lower it",
                    " in order to run analysis.")
            }
    }
    
    #adds the / Total#
    for (i in 1:length(rat_ttls)){
        incl_tbl[ names(rat_ttls)[i]] <- 
            paste(incl_tbl[names(rat_ttls)[i]], "/", 
                  rat_ttls[i], sep = "")
    }
    
    #prints the table of those included
    cat("Sample count: Include/Total\n")
    print(incl_tbl)
    #displays if flagged immediately
    if (low_samp_flag){
        warning(#"Sample Size Low:\n", 
            "One or more entries in above table have less included rats ",
            "than min_rat of ", min_rat, ".\n", immediate. = TRUE)
    }
    
    #formats the table so that information is easier to gleam from it
    prin_tbl <- as.data.frame(t_tbl %>%
                                  arrange(grepl(comp_vect[2], !!as.symbol(comp_vect[1])))) 
    #adds the rownames back
    for (i in 1:nrow(prin_tbl)){
        #adds more if they're both within else errors happen
        
        rownames(prin_tbl)[i] <- paste(prin_tbl[i, "ratID"], 
                                       prin_tbl[i, comp_vect[1]],
                                       sep = ":")
    }
    
    #prints the table
    cat("\n")
    cat("Individual Sample Nuclei Counts\n")
    print(prin_tbl)
    cat("\n")
    
    #rats that won't be in sample
    exclusion <- rownames(t_tbl[which(!t_tbl[,"Included"]),])
    
    #prints different message depending on # rats excluded
    if (length(exclusion) != 0){
        if (length(exclusion) > 1){
            cat(cat(exclusion, 
                    sep = ", "), 
                "to be excluded as counts < ", min_cell, "\n\n") 
        }
        else {
            cat(exclusion, " to be excluded as counts < ", min_cell, "\n\n")
        }
    } else {
        cat("No rats excluded\n\n")
    }
    
    #gets rid of the excluded rats
    temp <- t_tbl[which(t_tbl[,"Included"]),]
    
    #Checks to see if too many were excluded by the min_cell
    for (cmp in comp_vect[2:3]){
        if (!(cmp %in% temp[, comp_vect[1]])){
            stop("Exclusion excluded all rats with characteristic, ", 
                 cmp, ", for factor: ", comp_vect[1], "\n",
                 "  Please adjust min_cell to exclude less rats")
        }
    }
        
    #takes those rats out
    Idents(sub_obj) <- "ratID"
    sub_obj <- subset(sub_obj, idents = unique(temp$ratID))
    
    #just simply makes the table 
    counts <- to_pseudobulk(
        sub_obj, #The source of what we're generating a count
        replicate_col = "ratID",
        cell_type_col = "cluster_name",
        label_col = comp_vect[1]
    )[[1]][, rownames(temp)]
    
    #makes formula
    form_fact <- as.formula(paste("~", comp_vect[1], sep = " "))
    
    #turns matrix into DESeq object
    clust_tbl <- DESeqDataSetFromMatrix(counts, #Data Table
                                        colData = temp,  #metaData
                                        design = form_fact)
    
    # runs DESseq2
    clust_tbl <- DESeq(clust_tbl)
    
    # now to add option to save the DESeq object
    if (keep_dds == TRUE){
        #saveRDS(clust_tbl, file = paste("DESeq2_", cluster, ".rds", sep = ""))
      #assign(dds, dds, envir = .GlobalEnv)
    }
    
    #populates results
  
results = results(clust_tbl, 
                  contrast = comp_vect,
                  alpha = 0.05,
                  cooksCutoff = FALSE,
                  independentFiltering = FALSE)

    #turns into tibble
 results <- results %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
    
    #returns
    #return(results)
 return(list(deseq_results = results, clust_tbl = clust_tbl))
}
