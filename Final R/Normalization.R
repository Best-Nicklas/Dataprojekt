#' Description - Helper function to normalization of the constants.
#' 
#' @param MAF List of Minor Ale Frequency.
#' @param causal List of causal SNP´s where causal SNP´s has the value 1.
#' @return A list with the normalized values of mu and sigma.
#' @export

calc_normalization_consts <- function(MAF, causal) {
  mu <- 2 * MAF * causal
  sigma <- sqrt((2 * MAF * causal)*(1 - MAF * causal))
  sigma[sigma == 0] <- 1
  
  return(list(mu = mu, sigma = sigma))
}
