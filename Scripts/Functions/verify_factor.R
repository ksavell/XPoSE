#' Verify Factor
#' 
#' General function that will take in a vector of user input and correct it using
#' a vector of the wanted factor's names. 
#'
#' @param object SE
#' @param vect the vector being verified, should have at least 2 elems
#' @param factor the factor to verify with. Within our object this will be 
#'               'group' or 'sex'
#'
#' @return the user's (corrected) choice vector
#' @export
verify_factor <- function(object, vect, factor){
  #makes vector of valid options, normal and lowercase
  valids <- unique(object[[factor]])[, 1]
  valids_low <- tolower(valids)
  
  #There is a case where factor is the first elem so this ensures output to user
  #is consistent
  mod <- 0 
  
  #Loops vector
  for (i in 1:length(vect)){
    #Checks if first elem is the factor
    if (i == 1){
      if (vect[i] != factor){
        #Loops until valid input
        while (!(vect[i] %in% valids) & vect[i] != "d"){
          #checks for case, converts to proper if not
          if (tolower(vect[i]) %in% valids_low){
            vect[i] <- valids[match(tolower(vect[i]), valids_low)]
          }else {
            #prompts user about 
            cat("Invalid input for Comparison", i, 
                ". Please check for typo.\n",
                "Valid inputs include:\n",
                "If you would like to delete this option, please enter 'd'. \n",
                sep = "")
            sapply(valids, print, quote = FALSE)
            vect[i] = readline(paste("Comparsion", i, 
                                     ": ", sep = ""))
          }
        }
        #If the factor is the first index, i needs to be 1
      } else {
        mod <- 1
      }
    }else {
      #Loops until valid input
      while (!(vect[i] %in% valids) & vect[i] != "d"){
        #checks for case 
        if (tolower(vect[i]) %in% valids_low){
          vect[i] <- valids[match(tolower(vect[i]), valids_low)]
        }else {
          cat("Invalid input for Comparison", i - mod, 
              ". Please check for typo.\n",
              "Valid inputs include:\n",
              "If you would like to delete this option, please enter 'd'. \n",
              sep = "")
          sapply(valids, print, quote = FALSE)
          vect[i] = readline(paste("Comparsion", i - mod, ": ", sep = ""))
        }
      }
    }
  }
  
  #Check to see if we need to rerun the function
  rerun <- FALSE
  
  #dupe check
  if (length(vect) > 1){
    for (i in (1 + mod):(length(vect) - 1)){
      for(j in (2 + mod):(length(vect))){
        #Prevents issue where i and j can both be 2
        if (i != j){
          while (vect[i] == vect[j]){
            rerun <- TRUE
            cat("Dupe found at Comparsion", j - mod,
                ". Please use a different option within your set.\n",
                "Valid inputs include:\n", 
                "If you would like to delete this option, please enter 'd'. \n",
                sep = "")
            sapply(valids, print, quote = FALSE)
            cat("Please don't put in the same thing in...")
            vect[i + 1] = readline(paste("Comparsion", j - mod, ": ", sep = ""))
          }
        }
      }
    }
    #Cleans vect of d's
    vect <- vect[vect != "d"]
  }
  
  while (length(vect) < 2){
    rerun <- TRUE
    cat("Additional comparison needed.",
        " Please enter one of the below options:\n")
    sapply(valids, print, quote = FALSE)
    vect <- append(vect, readline("New Comparsion: "))
  }
  
  #Validates possible new changes
  if (rerun){
    vect <- verify_factor(object, vect, factor)
  }
  
  #returns vect with changes (if any)
  return(vect)
}
