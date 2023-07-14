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

        #prevents ratID use
        tbl_lst[["hidden"]] <- tbl_lst[["hidden"]][tbl_lst[["hidden"]] != "ratID"]

        #Users can append to hidden. This code prevents more than 2 factors
        #from existing in hidden when it is used. 
        if (length(tbl_lst[["hidden"]]) > 2){
          tbl_lst[["hidden"]] <- tbl_lst[["hidden"]][1:2]
        }
      
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
