#' Perform genome wide association study (GWAS)
#' 
#' This function performs a GWAS on given genotype and phenotype data to find out how statistically associated SNPs are with a target disease.
#' 
#' @param rds.obj A list of length 3 containing an FBM.code256 named genotypes and accompanying FAM and MAP tibbles.
#' @param y Vector of regressands to regress on (typically case-control status of the genotypes in rds.obj)
#' @param include Vector of rows to use in regression. Used with cross-validation. Default uses all rows.
#' @return A data.frame with slopes of each regression, standard errors of each slope, t-scores associated with each slope and P-values of each slope.
#' @export 
#' 
#' 
GWAS <- function(rds.obj, y, include = bigparallelr::rows_along(rds.obj$genotypes)) {
  FBM <- rds.obj$genotypes
  
  #Uses function from bigsnpr package to do regression on FBM
  regr <- bigstatsr::big_univLinReg(FBM, y[include], ind.train = include)
  
  #Adds column with P-values
  regr$p.value <- predict(regr, log10 = FALSE)
  
  return(data.frame(regr))
}
