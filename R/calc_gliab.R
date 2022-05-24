#' title
#' 
#' Helper function to calculate genetic liabilities.
#' @param obj placeholder
#' @param beta placeholder
#' @param mu placeholder
#' @param sigma placeholder
#' @return placeholder
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