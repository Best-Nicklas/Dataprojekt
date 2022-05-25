#' title Calculation for genetic liability
#' 
#' This is a helper function to used to calculate genetic liabilities.
#' 
#' @param obj A matrix where colums are genotypes and rows are amout of persons.
#' @param beta List of beta values.
#' @param mu List of mu vaules.
#' @param sigma List of sigma vaules.
#' @return The function returns a matrix with one genetic liability for each row
#' @keywords internal
#' @export

calc_gliab <- function(obj, beta, mu, sigma) {
  g_liab <- sweep(sweep(obj, 
                        MARGIN = 2, 
                        STATS = mu, 
                        FUN = "-"), 
                  MARGIN = 2, 
                  STATS = sigma, 
                  FUN = "/") %*% beta
  
  return(g_liab)
}