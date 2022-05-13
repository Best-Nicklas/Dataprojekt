#' Description - Helper function to calculate a covariance matrix
#' 
#' @param h2
#' @param sib
#' @return A covariance matrix, calculated by the value of h2 and 
#' the amount of siblings.
#' @example covmatrix(0.5, 2)
#' 

covmatrix <- function(h2, n_sib = 0) {
  cov <- matrix(h2/2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3,4] <- cov[4,3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  
  return(cov)
}