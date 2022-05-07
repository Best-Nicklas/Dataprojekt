#' Description - Helper function to calculate siblings liability 
#' 
#' @param n_sibs
#' @param p1
#' @param p2
#' @param disease
#' @return 
#' @example 
#' 

sibs_liab_calc <- function(n_sibs, p1, p2, disease) {
  if (n_sibs == 0) {
    return(NULL)
  }
  
  MAF <- disease$MAF
  beta <- disease$BETA
  causal <- disease$CAUSAL
  
  #Generates children from p1 and p2
  temp <- truer_gen(p1, p2, NULL)
  
  #Calculates mu and sigma used for normalization
  mu <- 2 * MAF * causal
  sigma <- sqrt((2 * MAF * causal)*(1 - MAF * causal))
  sigma[sigma == 0] <- 1
  
  c <- sum(causal)
  
  #Calculates the genetic liabilities for sibs
  liab_g <- c(sweep(sweep(temp, 
                          MARGIN = 2, 
                          STATS = mu, 
                          FUN = "-"), 
                    MARGIN = 2, 
                    STATS = sigma, 
                    FUN = "/") %*% beta)
  return(liab_g)  
}

