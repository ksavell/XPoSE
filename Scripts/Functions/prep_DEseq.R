#' Get Pie Pieces
#' 
#' Returns chosen cluster's count table
#'
#' @param data_tab Aggregate table for the given object, returned by 
#'                 'get_cnt_tbls'
#' @param data_lst List of object's split cluster objects, returned by 
#'                 'make_cluster_list'
#' #@param clust_vect Vector of used clusters, returned by 'Kowalski'
#' @param object SE
#' @param dataname name of object, used a lot for once, if you put in the wrong
#'                 name, this function breaks
#'
#' @return the user's chosen count table
#' @export
#'
#' @examples
#' #Get count table for a cluster in the glut object
#' my_cnt_tbl <- prep_DEseq(glut_tab, glut, "Glut")
prep_DEseq <- function(data_tab, object, dataname){
  # prompt user for table
  cat("What cluster do you want to see a count table for?", 
      "Press 'v' to view options.", sep = "\n")
  piece = readline()
  
  # vectors to store cluster names
  data_lst <- c()
  clust_vect <- c()
  
  # populates data_lst 
  for (name in names(unique(data_tab))){
    if (str_detect(name, dataname)){
      data_lst <- append(data_lst, name)
    }
  }
  
  # splits data_lst into a list 
  temp <- sapply(data_lst, strsplit, dataname)
  
  # populates clust_vect using numbers from data_lst
  for (i in 1:(length(temp))){
    clust_vect <- append(clust_vect, temp[[i]][2]) 
  }  
  
  # removes possible NAs from clust_vect
  clust_vect <- clust_vect[!is.na(clust_vect)]  
  
  # Verifies 
  while (!((piece %in% data_lst) | (piece %in% clust_vect)) 
         & piece != "all"){
    if (piece != "v"){
      cat("Invalid input. Please enter a valid cluster or 'v'. \n") 
    }
    # Shows options
    cat("Your options are:\n")
    sapply(data_lst, print, quote = FALSE)
    sapply(clust_vect, print, quote = FALSE)
    print("all", quote = FALSE)
    print("v", quote = FALSE)
    # reruns prev line
    cat("What cluster do you want to see a count table for?", 
        "Press 'v' to view options.", sep = "\n")
    piece = readline()
  }
  
  # if a single num (like a 4), turns to something we can check for in the data Table
  if (piece %in% clust_vect){
    piece <- paste(dataname, piece, sep = "")
  }
  # same as prev but with all
  if (piece == "all"){
    piece <- paste("all", dataname, sep = "")
  }
  
  # Finds section
  tbl_sect <- match(piece, names(data_tab))
  
  # gets piece of the pie (table)
  tbl_piece <- data_tab[1:nrow(data_tab), (tbl_sect + 1):(tbl_sect + length(unique(object@meta.data[["ratID"]])))]
  return(tbl_piece)
}
