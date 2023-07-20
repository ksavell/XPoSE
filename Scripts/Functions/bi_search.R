#' Binary Search
#' 
#' Recursively finds target elem's index within the container using a binary 
#' search algorithm
#' If we release this in a package, I will need to verify it's an atomic vector
#'
#' @param cntnr an atomic vector containing the target
#' @param target elem of cntnr you're trying to find index of
#' @param left first index of wanted range (default is 1)
#' @param right last index of wanted range (default is length(cntnr))
#'
#' @return the index of the target
#' @export
#'
#' @examples
#' #Find index of "O" in LETTERS
#' bi_search(LETTERS, "O")
#' #Find index of "O" in a specific range (27-53)
#' bi_search(c(LETTERS, LETTERS, LETTERS), "O", 27, 53)
bi_search <- function(cntnr, target, left = 1, right = length(cntnr)){
  
  if (right >= 2 & left <= right){
    #finds middle
    mid <- as.integer((right + left)/2)
    
    #target found
    if (cntnr[mid] == target){
      return(mid)
    
    #target greater
    } else if (cntnr[mid] < target){
      return(bi_search(cntnr, target, mid + 1, right))
      
    #target lesser
    } else {
      return(bi_search(cntnr, target, left, mid - 1))
    }
    
  #object not found
  }else {
    return(NA)
  }
}
