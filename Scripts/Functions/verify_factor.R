#' Verify Factor
#' 
#' General function that will take in a vector of user input and correct it using
#' a vector of the wanted factor's names. 
#'
#' @param object SE
#' @param vect the vector being verified
#' @param factor the factor to verify with. Within our object this will be 
#'               'group' or 'sex'
#'
#' @return the user's (corrected) choice vector
#' @export
verify_factor <- function(object, vect, factor){
  
  #makes vector of valid options
  valids <- unique(object[[factor]])[,1]
  
  #Loops vector
  for (i in 1:length(vect)){
    #There is a case where factor is the first elem so this ensure output to user
    #is consistent
    mod <- 0 
    if (i == 1){
      if (vect[i] != factor){
        #Loops until valid input
        while (!(vect[i] %in% valids)){
          cat("Invalid input for Comparison", i, ". Please check for typo.\n",
              "Valid inputs include:\n", sep = "")
          sapply(valids, print, quote = FALSE)
          vect[i] = readline(paste("Comparsion", i, ": ", sep = ""))
        }
        #If the factor is the first index, i needs to be 1
      } else {
        mod <- 1
      }
    }else {
      #Loops until valid input
      while (!(vect[i] %in% valids)){
        cat("Invalid input for comparison", (i + mod), ". Please check for typo.\n",
            "Valid inputs include:\n", sep = "")
        sapply(valids, print, quote = FALSE)
        vect[i] = readline(paste("comparison", (i + mod), ": ", sep = ""))
      }
    }
  }
  
  #returns vect with changes (if any)
  return(vect)
}