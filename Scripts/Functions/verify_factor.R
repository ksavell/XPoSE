#' Verify Factor
#' 
#' General function that will take in a vector of user input and correct it using
#' a vector of the wanted factor's names. 
#'
#' @param dataset SE
#' @param vect the vector being verified
#' @param factor the factor to verify with. Within our dataset this will be 
#'               'group' or 'sex'
#'
#' @return the user's (corrected) choice vector
#' @export
verify_factor <- function(dataset, vect, factor){
  
  #makes vector of valid options
  valids <- unique(dataset[[factor]])[,1]
  
  #Loops vector
  for (i in 1:length(vect)){
    #There is a case where group is the first elem so this ensure output to user
    #is consistent
    mod <- 0 
    if (i == 1){
      if (vect[i] != factor){
        #Loops until valid input
        while (!(vect[i] %in% valids)){
          cat("Invalid input for Group", i, ". Please check for typo.\n",
              "Valid inputs include:\n", sep = "")
          sapply(valids, print, quote = FALSE)
          vect[i] = readline(paste("Group", i, ": ", sep = ""))
        }
        #If group = 1, i needs subtraction
      } else {
        mod <- 1
      }
    }else {
      #Loops until valid input
      while (!(vect[i] %in% valids)){
        cat("Invalid input for Group", (i + mod), ". Please check for typo.\n",
            "Valid inputs include:\n", sep = "")
        sapply(valids, print, quote = FALSE)
        vect[i] = readline(paste("Group", (i + mod), ": ", sep = ""))
      }
    }
  }
  
  #returns vect with changes (if any)
  return(vect)
}