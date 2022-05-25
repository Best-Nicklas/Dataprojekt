#' Truncated
#' 
#' Helper function to calculate threshold for truncated distribution.
#' 
#' @param n Integer
#' @param range placeholder
#' @param mu Mu
#' @param sigma Sigma
#' @return The threshold value used in a truncated distribution
#' @export
#' 

rnorm_trunc <- function(n, range, mu, sigma) {
  
  lower <- pnorm(min(range), mu, sigma)
  upper <- pnorm(max(range), mu, sigma)
  
  u <- runif(n, lower, upper)
  
  return(qnorm(u, mu, sigma))
}