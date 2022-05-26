#' Calculation of genetic liabilities
#' 
#' This is an internal helper function to used to calculate genetic liabilities.
#' 
#' @param obj A matrix consisting of 0s, 1s or 2s, where each row represents a genotype.
#' @param beta List of beta values.
#' @param mu List of mu values to normalize on.
#' @param sigma List of sigma values to normalize on.
#' @return The function returns a vector of genetic liabilities for each genotype.
#' @keywords internal
#' @export

calc_gliab <- function(obj, beta, mu, sigma) {
  if (ncol(obj) != length(beta)) stop("Number of columns in obj must be equal to length of beta")
  if (ncol(obj) != length(mu)) stop("Number of columns in obj must be equal to length of mu")
  if (ncol(obj) != length(sigma)) stop("Number of columns in obj must be equal to length of sigma")
  
  # Uses sweep to normalize obj using the normalization constants mu and sigma, then does 
  # matrix multiplication with beta.
  g_liab <- sweep(sweep(obj, 
                        MARGIN = 2, 
                        STATS = mu, 
                        FUN = "-"), 
                  MARGIN = 2, 
                  STATS = sigma, 
                  FUN = "/") %*% beta
  
  return(g_liab)
}

