#' Make Mega Table
#'
#' @param object SE -- drakeeeee
#'
#' @return a 'MEGA' table of all of the chosen cluster's count tables
#' @export
#'
#' @examples
#' #makes an aggragate table combining all the count tables of user chosen 
#' #clusters
#' my_glut_cnts <- mk_mega_tbl(glut)
mk_mega_tbl <- function(object){
  
  #prompts user for what they want their data to be named
  cat("What would you like the name for your object to be?",
      "Good examples include:",
      "Glut",
      "GABA", sep = "\n")
  dataname = readline()
  
  #confirms the name
  print(paste("Are you sure you want your object to be named '", dataname, "'",
              "? (y/n)", sep = ""), quote = FALSE)
  confirm = readline()
  
  #gives user chance to change
  while (confirm != "y"){
    #invalid input
    if (confirm != "n"){
      confirm = readline("Invalid input. Please enter a 'y' or 'n'. ")
    }else {
      #prompts user again
      cat("What would you like the name for your object to be?",
          "Good examples include:",
          "Glut",
          "GABA", sep = "\n")
      dataname = readline()
      print(paste("Are you sure you want your object to be named '", dataname, "'",
                  "? (y/n)", sep = ""))
      confirm = readline()
    }
  }
  
  #Makes our data containers
  # clust_nams <- c()
  # clust_lst <- list()
  data_tabs <- matrix()
  
  #gets and stores the cluster variables that user defines
  #clust_nams <- Kowalski(object, dataname)
  
  #creates and stores seurot objects from split clusters
  #clust_lst <- make_clust_list(clust_nams, object, dataname)
  
  #makes the covetted MEGA table
  #data_tabs <- get_cnt_tbls(clust_lst, clust_nams, object, dataname)
  data_tabs <- get_cnt_tbls(make_clust_list(clust_nams, object, dataname),
                            Kowalski(object, dataname),
                            object,
                            dataname)
  
  #returns the MEGA table
  return(data_tabs)
}
