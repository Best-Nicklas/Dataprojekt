#' Helper function to a covariance matrix
#'
#' @param h2 Heritability parameter
#' @param sib Number of siblings 
#' @return A covariance matrix for the liabilities of a family with /code{sib}
#' number of siblings
#' @examples   
#' covmatrix(0.8, 2)

covmatrix <- function(h2, sib = 0) {
  cov_m <- matrix(c(h2, h2, rep(0.5 * h2, sib + 2),
                    h2, 1, rep(0.5 * h2, sib + 2),
                    0.5 * h2, 0.5 * h2, 1, 0, rep(0.5 * h2, sib),
                    0.5 * h2, 0.5 * h2, 0, 1, rep(0.5 * h2, sib)),
                  nrow = 4, ncol = 4 + sib, byrow = T)
  
  if (sib != 0) {
    for (i in 1:sib) {
      cov_m <- rbind(s, c(rep(0.5 * h2, 3 + i), 1, rep(0.5 * h2, sib - i)))
    }
  }
  
  return(cov_m)
}