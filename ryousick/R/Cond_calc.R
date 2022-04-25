#' Helper function to calculate mu and sigma
#'
#' @param i A integer indicating which index of a liability we want to
#' calculate the distribution for.
#' @param covmatrix A covariance matrix of a multivariate normal distribution.
#' @return The a list with the newly calculated mu and sigma
#' @examples   
#' cond_calc(1, covmatrix)


cond_calc <- function(i, covmatrix) {
  s11 <- covmatrix[i, i]
  s12 <- covmatrix[i, -i]
  s21 <- covmatrix[-i, i]
  s22 <- covmatrix[-i, -i]
  
  new_mu <- s12 %*% solve(s22) 
  new_sigma <- s11 - (s12 %*% solve(s22) %*% s21)
  
  return(c(new_mu, new_sigma))
}