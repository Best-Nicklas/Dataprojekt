#' Description - Helper function to calculate a covariance matrix
#' 
#' @param h2
#' @param sib
#' @return A covariance matrix, calculated by the value of h2 and 
#' the amount of siblings.
#' @example covmatrix(0.5, 2)
#' 

covmatrix <- function(h2, sib = 0) {
  cov_m <- matrix(c(h2,       h2, rep(0.5 * h2, sib + 2),
                    h2,        1, rep(0.5 * h2, sib + 2),
                    0.5 * h2, 0.5 * h2, 1, 0, rep(0.5 * h2, sib),
                    0.5 * h2, 0.5 * h2, 0, 1, rep(0.5 * h2, sib),
                    unlist(purrr::map(1:sib, .f = ~ 
                                        {c(rep(0.5 * h2, 3 + .x), 1, rep(0.5 * h2, sib - .x))}))),
                  nrow = 4 + sib, ncol = 4 + sib, byrow = T)
  
  return(cov_m)
}

# det her er en meget bøvlet måde at lave denne matrix på.
# udnyt i ved hvad de fleste indgange er, og så tilpas resten.
get_cov = function(h2, n_sib = 0) {
  cov <- matrix(h2/2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3,4] <- cov[4,3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  cov
}
