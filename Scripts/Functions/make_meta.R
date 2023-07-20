#' Make Metadata
#' 
#' Swipes the rat:group pairs from the big table in order to make the metadata
#'
#' @param data_tab 
#'
#' @return returns metadata as data frame
#' @export
#'
#' @examples
#' #Make metadata for either object (glut_tabs is properly populated)
#' metaData <- make_meta(glut_tabs)
make_meta <- function(data_tab){
  #This *should* work for another experiment as long as the big table has at 
  #least one cluster in it.
  ratG <- names(data_tab[2:9])
  
  #splits into list and makes storage for just the groups
  temp <- sapply(ratG, strsplit, ":")
  justG <- c()
  
  #Fills justG w/ the groups in the split
  for (i in 1:length(temp)){
    justG <- append(justG, temp[[i]][2]) 
  }
  
  #Makes data frame and labels it 
  met_fram <- data.frame(justG)
  rownames(met_fram) <- ratG
  colnames(met_fram) <- "group"  #<- Will likely need to be changed to a vector if sex also wanted
  
  #checks that the metadata is set up correctly
  cat("There should be a 'TRUE' under this line:\n")
  cat(as.character(all(rownames(met_fram) == colnames(counts)), "\n", sep = ""))
  cat("\n")
  
  #returns frame
  return(met_fram)
}
