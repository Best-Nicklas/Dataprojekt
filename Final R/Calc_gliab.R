#' Description - Helper function to calculate genetic liabilities.
#' @param obj 
#' @param beta
#' @param mu
#' @param sigma
#' @return 
#' @example 
#' 

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