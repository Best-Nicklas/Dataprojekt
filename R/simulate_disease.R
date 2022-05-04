#' Description - 
#' 
#' @param n_snp
#' @param prevalens
#' @param h2
#' @param causal
#' @param maf
#' @param maf_low
#' @param maf_high
#' @param seed
#' @return 
#' @example 
#' 


simulate_disease <- function(n_snp, prevalens, h2, causal = NULL, maf = NULL, maf_low = 0.01, maf_high = 0.49, seed = NULL){
  if (!is.null(seed)) {set.seed(seed)}
  
  if (is.null(causal)) {
    causal <- numeric(n_snp)
    causal[sample(1:n_snp, max(1, round(n_snp / 10)))] <- 1
  }
  
  if (is.null(maf))  {
    maf <- runif(n_snp, maf_low, maf_high)
    
  }
  
  c <- sum(causal)
  beta <- ifelse(causal == 1, rnorm(n_snp, 0, sqrt(h2 / c)), 0)
  
  return(list(MAF = maf, BETA = beta, CAUSAL = causal, H2 = h2, PREVALENS = prevalens, N_SNP = n_snp))
}
