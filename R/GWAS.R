#' Description - GWAS
#' 
#' @param person
#' @param y
#' @return 
#' @example 
#' 


GWAS <- function(person, y) {
  FBM <- person$genotypes
  
  regr <- big_univLinReg(FBM, y)
  regr$p.value <- predict(regr, log10 = FALSE)
  
  return(data.frame(pval = regr$p.value))
}
