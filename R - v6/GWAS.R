#' Description - GWAS
#' 
#' @param person
#' @param y
#' @return 
#' @example 
#' 


GWAS <- function(person, y) {
  
  FBM <- person$genotypes
  
  #uses function from bigSNPr package to do regression on FBM
  regr <- big_univLinReg(FBM, y)
  #adds column with pvalues
  regr$p.value <- predict(regr, log10 = FALSE)
  
  return(data.frame(regr))
}
