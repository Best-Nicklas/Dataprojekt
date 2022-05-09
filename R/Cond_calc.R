#' Description - Helper function to calculated the conditional mu and sigma
#' 
#' @param i
#' @param covmatrix
#' @return The conditional values of mu and sigma
#' @example 
#' 


cond_calc <- function(i, covmatrix) {
  s11 <- covmatrix[i, i]
  s12 <- covmatrix[i, -i]
  s21 <- covmatrix[-i, i]
  s22 <- covmatrix[-i, -i]
  
  new_mu <- s12 %*% solve(s22) 
  new_sigma <- s11 - (s12 %*% solve(s22) %*% s21)
  
  return(c(new_mu, new_sigma)) # jeg er IKKE fan af at lave dette trick her. 
}
#i den returnerende vektor, kan i NEMT løbe ind i problemer med at uforudsete ting sker. lav i stedet en liste med mu og sigma som navngivet indgange.

