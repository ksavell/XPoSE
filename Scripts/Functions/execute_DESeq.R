#' Execute DESeq
#'
#'Converts a list of count tables into a list of result tibbles
#'
#' @param tbl_lst a list of count tables generated from cluster data
#' @param object the Seurot object to do analysis on
#' @param factor the factor to verify with. Within our object this will be 
#'               'group' or 'sex' 
#' @param comp_vect Aspects within the factor you would like to compare. First 
#'                  elem will become your Experimental Variable, the second will
#'                  become your control
#'
#' @return a list of results tables generated from the count tables
#' @export
#'
#' @examples
execute_DESeq <- function(tbl_lst, object, factor, comp_vect){
  
  #splits names into list of the ratID's and groups
  temp <- sapply(names(tbl_lst[[1]]), strsplit, ":")
  justF <- c()
  
  #Fills justF w/ only the comparisons in the split
  for (i in 1:length(temp)){
      justF <- append(justF, temp[[i]][2]) 
  }
  
  #Makes data frame and labels it 
  met_fram <- data.frame(justF)
  rownames(met_fram) <- names(tbl_lst[[1]])
  colnames(met_fram) <- factor 
  
  #checks that the metadata is set up correctly
  cat("There should be at least one 'TRUE' under this line:\n")
  cat(as.character(unique(met_fram[,1]) %in% unique(object@meta.data[factor][, 1])), "\n")
  cat("\n")
  
  #makes necessary materials from params
  contrast_vect <- c(factor, comp_vect)
  form_fact <- as.formula(paste("~", factor, sep = ""))
  
  #our base list to hold tibbles, will be returned at function's end
  agg_res <- list()
  
  #Makes result table for each count table in list (excludes hidden)
  for (i in 1:length(tbl_lst)){
      #Makes object
      cluster <- DESeqDataSetFromMatrix(tbl_lst[[i]],  #Data Table
                                        colData = met_fram,  #metaData
                                        design = form_fact)
      
      #Filters out insignificant parts
      cluster <- cluster[ rowSums(counts(cluster)) > 5, ]
      cluster <- DESeq(cluster)
      objectname <- 
      assign(paste0(names(tbl_lst)[i],"_dds"), cluster, envir = .GlobalEnv)
      
      #Makes results obj
      results <- results(cluster,
                         contrast = contrast_vect,
                         alpha = 0.05,
                         cooksCutoff = FALSE,
                         independentFiltering = FALSE)
      
      #Turns results into a Tibble
      results_tib <- results %>%
                     data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     as_tibble()
      
      #appends to list of tibbles
      agg_res <- list.append(agg_res, clust_nam = results_tib)
      names(agg_res)[i] <- names(tbl_lst)[i]
  }
  
  #returns results
  return(agg_res)
}
