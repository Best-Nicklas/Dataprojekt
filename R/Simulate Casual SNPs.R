#' Description 
#'
#' @param n 
#' @param c 
#' @return 
#' @examples   
#' Simulate_Causal_SNPs(0.8, 2, 777)

Simulate_Causal_SNPs <- function(n, c, seed = NULL){
  if (!is.null(seed)) {set.seed(seed)}
  Causal_SNP <- sample(c(1:n),c)
  return(Causal_SNP)
}
