#' Description - Helper function to calculated the interval for a truncated distribution
#' 
#' @param n
#' @param range
#' @param mu
#' @param sigma
#' @return A interval which can be use to pull from disstribution truncated.
#' @example 
#' 

rnorm_trunc <- function(n, range, mu, sigma) {
  
  lower <- pnorm(min(range), mu, sigma)
  upper <- pnorm(max(range), mu, sigma)
  
  u <- runif(n, lower, upper)
  
  return(qnorm(u, mu, sigma))
}
