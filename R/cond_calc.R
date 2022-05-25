#' Conditional values of mu and sigma
#' 
#' This function is an internal helper function used to calculate the 
#' conditional mu and sigma of a co-variance for a multivariate normal distribution.
#'  
#' @param i parameter index of interest.
#' @param covmatrix the co-variance matrix. 
#' @return The conditional values of mu and sigma.
#' @keywords internal
#' @export

cond_calc <- function(i, covmatrix) {
  s11 <- covmatrix[i, i]
  s12 <- covmatrix[i, -i]
  s21 <- covmatrix[-i, i]
  s22 <- covmatrix[-i, -i]
  
  new_mu <- s12 %*% solve(s22) 
  new_sigma <- s11 - (s12 %*% solve(s22) %*% s21)
  return(list(mu = new_mu, sigma = new_sigma))
}
