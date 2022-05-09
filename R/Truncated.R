#' Description - Helper function to get values from truncated normal distribution
#' 
#' @param n
#' @param range
#' @param mu
#' @param sigma
#' @return Values from a normal distribution truncated in the given range.
#' @example 
#' 

rnorm_trunc <- function(n, range, mu, sigma) {
  
  #hvis i bare bruger 0 til 1, som jeres normale skala, så skal i ikke skifte mellem reelle tal og (0,1) og så tilbage til reelle tal igen. (for grænserne)
  lower <- pnorm(min(range), mu, sigma)
  upper <- pnorm(max(range), mu, sigma)
  
  u <- runif(n, lower, upper)
  
  return(qnorm(u, mu, sigma))
}
