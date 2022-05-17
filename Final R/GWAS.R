#' Description - GWAS
#' 
#' @param person
#' @param y
#' @param include
#' @return 
#' @example 
#' 


GWAS <- function(person, y, include = rows_along(person$genotypes)) {
  FBM <- person$genotypes
  
  #uses function from bigSNPr package to do regression on FBM
  regr <- big_univLinReg(FBM, y[include], ind.train = include)
  #adds column with pvalues
  regr$p.value <- predict(regr, log10 = FALSE)
  
  return(data.frame(regr))
}
