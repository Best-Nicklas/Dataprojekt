#' Description - Function to calculate means of genetic liability 
#' from a configarations of famalies
#' 
#' @param config A configarations 
#' @param n Iterations after burn_in
#' @param start Iterations until mean and sigma convegates 
#' @param cov_mat Covariance matrix
#' @return The mean of genetic liability from the given configaration 
#' @example 
#'

gibbs_sampler <- function(config, n, start, cov_mat) {
  l_n <- dim(cov_mat)[1] #number of l's
  threshold <- qnorm(0.95)
  
  gen_liabs <- numeric(n)
  liabs_current <- rep(10,l_n) #initializing l's
  
  #pre-calculations for each l
  means <- matrix(ncol = l_n - 1, nrow = l_n)
  sigmas <- vector()
  for (p in 1:l_n) {
    temp <- cond_calc(p, cov_mat) 
    means[p, ] <- temp[1:(l_n - 1)]
    sigmas[p] <- temp[l_n]
  }
  
  # iterations
  for (i in 1:n) {
    # iteration for each parameter
    for (p in 1:l_n) {
      new_mean <- means[p, ] %*% liabs_current[-p]
      
      #For genetic liability - no truncation 
      if (p == 1) {
        liabs_current[p] <- rnorm(1, new_mean, sqrt(sigmas[p]))
      }
      
      #For liabilties when we have case (0)
      else if (config[p-1] == 0) {
        liabs_current[p] <- rnorm_trunc(1, c(-Inf, threshold), new_mean, sqrt(sigmas[p]))
        
      }
      #For liability when we dont have case (1)
      else {
        liabs_current[p] <- rnorm_trunc(1, c(threshold, Inf), new_mean, sqrt(sigmas[p]))
      }
    }
    gen_liabs[i] <- liabs_current[1]
  }
  return(mean(gen_liabs[start:n]))
}