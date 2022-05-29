#' Simulation of disease parameters
#' 
#' This function is used to simulate disease parameters. If not predetermined,
#' the function will randomly calculate causal SNPs, MAF vaules and beta values. 
#' 
#' @param n_snp Integer specifying amount of SNPs.
#' @param prevalence Integer specifying prevalence of disease in simulated population.
#' @param h2 Integer specifying heritability parameter.
#' @param causal Vector of predetermined causal SNPs (0 if not causal, 1 if causal). Leave empty to get random causal SNPs. 
#' @param causal_n Integer specifying amount of causal SNPs if no causal param given. Default value is 10 procent of total SNPs.
#' @param maf Vector of predetermined Minor Allele Frequencies.Leave empty to get random MAF.
#' @param maf_low Integer specifying lower bound for MAF if none given. 
#' @param maf_high Integer specifying upper bound for MAF if none given. 
#' @return A list with all the disease parameters.
#' @export
#' 


sim_disease <- function(n_snp, prevalence, h2, causal = NULL, causal_n = max(1,round(n_snp / 10)), maf = NULL, maf_low = 0.01, maf_high = 0.49){
  if (n_snp <= 0) stop("n_snp must be positive")
  if (causal_n < 0) stop("causal_n must be non-negative")
  if (prevalence <= 0 || prevalence >= 1) stop("prevalence must be between 0 and 1")
  if (h2 <= 0 || h2 >= 1) stop("h2 must be between 0 and 1")
  if (!is.null(causal) && length(causal) != n_snp) stop("length of causal must be equal to n_snp")
  if (!is.null(maf) && length(maf) != n_snp) stop("length of MAF must be equal to n_snp")
  if (maf_low < 0 || maf_low >= 1 || maf_low > maf_high) stop("maf_low must be between 0 and 1 and less than or equal maf_high")
  if (maf_high <= 0 || maf_low > 1 || maf_low > maf_high) stop("maf_high must be between 0 and 1 and greater than or equal maf_high")
    
  #Calculates a vector with causal SNPs at random positions if a causal vector not given
  if (is.null(causal)) {
    causal <- numeric(n_snp)
    causal[sample(1:n_snp, max(1, causal_n))] <- 1
  }
  
  #Calculates MAF if not given
  if (is.null(maf))  {
    maf <- runif(n_snp, maf_low, maf_high)
    
  }
  #Calcualtes betas for causal SNPs
  c <- sum(causal)
  beta <- ifelse(causal == 1, rnorm(n_snp, 0, sqrt(h2 / c)), 0)
  
  return(list(MAF = maf, BETA = beta, CAUSAL = causal, H2 = h2, PREVALENCE = prevalence, N_SNP = n_snp))
}
