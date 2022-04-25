#' Description 
#'
#' @param n 
#' @param c 
#' @return 
#' @examples   
#' Simulate_Causal_SNPs(n = 1000, c = 100, seed = 42)

Simulate_Causal_SNPs <- function(n, c, seed = NULL){
  if (!is.null(seed)) {set.seed(seed)}
  causal_SNP_loc <- sample(c(1:n),c)
  causal_SNP <- numeric(n)
  causal_SNP[causal_SNP_loc] <- 1
  return(causal_SNP)
}