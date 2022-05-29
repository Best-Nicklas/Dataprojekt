#' Calculate co-variance matrix for liability modeling
#' 
#' This function creates the theoretical co-variance matrix for the 
#' multivariate normal distribution which is used to model the liabilities.
#' 
#' @param h2 Integer specifying heritability parameter (0 < h2 < 1).
#' @param n_sib Integer specifying number of siblings (default 0).
#' @return A co-variance matrix, created from the value of h2 and 
#' the number of sibs. 
#' 
#' @export
#' 
covmatrix <- function(h2, n_sib = 0) {
  if (n_sib < 0) stop("n_sib must be a positive integer")
  if (h2 >= 1) stop("h2 must be an integer between 0 and 1")
  if (h2 <= 0) stop("h2 must be an integer between 0 and 1")
  
  # the base structure
  cov <- matrix(h2/2, 4 + n_sib, 4 + n_sib)
  # diagnoal entries 1
  diag(cov) <- 1
  # positions where 0
  cov[3,4] <- cov[4,3] <- 0
  # positions where h2
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  
  return(cov)
}