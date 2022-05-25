#' Calculation of normalization constants
#' 
#' This an internal helper function used to calculate the theoretical value of the constants 
#' mu and sigma, which are used to normalize the data.
#'  
#' @param MAF Vector of Minor Allele Frequencies.
#' @param causal Vector of causal SNPs where causal SNPs has the value 1 and non causal SNPs has the value 0.
#' @return A list with that contains two vectors, one for theoretical values of mu and one for the theoretical values of 
#' sigma.
#' @keywords internal
#' @export

calc_normalization_consts <- function(MAF, causal) {
  mu <- 2 * MAF * causal
  sigma <- sqrt((2 * MAF * causal)*(1 - MAF * causal))
  sigma[sigma == 0] <- 1
  
  return(list(mu = mu, sigma = sigma))
}
