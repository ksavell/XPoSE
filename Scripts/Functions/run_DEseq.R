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
