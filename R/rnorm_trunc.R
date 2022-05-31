#' Sample from a truncated normal distribution
#' 
#' Helper function used to sample from a truncated normal distribution.
#' 
#' @param n Integer value specifying the number of samples.
#' @param range vector with two values specifying the truncated range (lower bound and upper bound). 
#' @param mu Integer value specifying the mean of the normal distribution.
#' @param sigma Integer value specifying the standard deviation of the normal distribution.
#' @return Sample values from the truncated normal distribution.
#' @export
#' 

rnorm_trunc <- function(n, range, mu, sigma) {
  if(n <= 0) stop("n must be positive")
  
  lower <- pnorm(min(range), mu, sigma)
  upper <- pnorm(max(range), mu, sigma)
  
  u <- runif(n, lower, upper)
  
  return(qnorm(u, mu, sigma))
}