

make_volcano <- function(res_tbl, object, hover = FALSE, up_color = "red",
                         down_color = "blue", export = "FALSE"){
  library(dplyr)
  
  #makes new metadata
  obj <- object
  obj$cluster_group <- paste(obj$seurat_clusters, obj$group)
  
  #
  #tibble being used 
  if (length(class(res_tbl)) > 1 & class(res_tbl) == "tbl_df"){
    #tibble
    #if str_dtect > 1
    
  
  } else {
    #list being used
    if (class(res_tbl) == "list"){
      #list
      if (length(res_tbl) > 1){
        #tell user to pick elem in the list 
        
        #use names and match, 
        #return(make_volcano(res_tbl[[match]], object))
        
      #single elem list   
      } else {
        #tbl <- res_tbl[[1]]
        return(make_volcano(res_tbl[[1]], object))
      }
      
    #data.frame being used
    } else if (class(res_tbl) == "data.frame"){
      #data frame
      
      #check for multiples using str_detect tech
      #return if length > 1
    
    #not supported class 
    } else {
      #error message
      stop("res_tbl, ", deparse(substitute(res_tbl)), ", is not a supported ",
           "class. Please use an object of class 'tibble', 'list', or ", 
           "'data.frame'")
    }
  }
}


