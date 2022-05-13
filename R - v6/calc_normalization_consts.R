#' Description - Helper function for calculating normalization constants 
#' 
#' @param n
#' @param range
#' @param mu
#' @param sigma
#' @return 
#' @example 
#' 


calc_normalization_consts <- function(MAF, causal) {
  mu <- 2 * MAF * causal
  sigma <- sqrt((2 * MAF * causal)*(1 - MAF * causal))
  sigma[sigma == 0] <- 1
  
  return(list(mu = mu, sigma = sigma))
}