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