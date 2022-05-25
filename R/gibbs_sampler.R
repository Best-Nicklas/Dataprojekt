#' Gibbs sampling
#' 
#' This function is used as a helper function for the LTFH function. It calculates the 
#' posterior mean genetic liability for and individual with a given 
#' configuration (case-control status of individual, 
#' parents and siblings) when a co-variance matrix modeling the
#' liabilities is given.  
#' 
#' @param config vector of configuration, with 1 when case and 0 when not case. First position specifies
#' case-control status of target individual, second position that of parent 1, third of parent 2, 
#' and the rest that of siblings.
#' @param burn_in An integer that specifies the burn-in period before sampling for the mean.
#' @param cov_mat Co-variance matrix modeling the liabilities.
#' @param prevalence Integer specifying the likelihood of having the disease in the population.
#' @return The mean posterior genetic liability of the configuration.
#' @export
#'

gibbs_sampler <- function(config, burn_in, cov_mat, prevalence) {
  l_n <- nrow(cov_mat) #number of l's
  threshold <- qnorm(prevalence, lower.tail = F) 
  
  gen_liabs <- numeric(burn_in + 10000)
  liabs_current <- rep(10,l_n) #initializing l's
  
  #pre-calculations for each l
  means <- matrix(ncol = l_n - 1, nrow = l_n)
  sigmas <- vector()
  for (p in 1:l_n) {
    temp <- cond_calc(p, cov_mat) 
    means[p, ] <- temp$mu
    sigmas[p] <- temp$sigma
  }
  
  SEM <- 1
  i <- 0
  
  # iterations
  while (SEM > 0.01) {
    i <- i + 1
    # iteration for each parameter
    for (p in 1:l_n) {
      new_mean <- means[p, ] %*% liabs_current[-p]
      
      #For genetic liability - no truncation 
      if (p == 1) {
        liabs_current[p] <- rnorm(1, new_mean, sqrt(sigmas[p]))
      }
      
      #For liabilties when we dont have case (0)
      else if (config[p-1] == 0) {
        liabs_current[p] <- rnorm_trunc(1, c(-Inf, threshold), new_mean, sqrt(sigmas[p]))
        
      }
      #For liability when we do have case (1)
      else {
        liabs_current[p] <- rnorm_trunc(1, c(threshold, Inf), new_mean, sqrt(sigmas[p]))
      }
    }
    gen_liabs[i] <- liabs_current[1]
    
    if (i > burn_in + 1000){
      SEM <- sd(gen_liabs[burn_in:i]) / sqrt(i - burn_in)
      
    }
  }
  return(mean(gen_liabs[burn_in:i]))
}

