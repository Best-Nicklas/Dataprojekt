#' Description
#'
#' @param n 
#' @param low_freq
#' @param high_freq
#' @param seed 
#' @return 
#' @examples   

SimulateMAF <- function(n, low_freq = 0.01, high_freq = 0.49, seed = NULL){
  if (!is.null(seed)) {set.seed(seed)}
  return(runif(n,low_freq, high_freq))
}