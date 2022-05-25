#' Perform GWAS with proxy information on family case-control status
#' 
#' This function uses proxy information on case-control status and given genotype data to find the likelihood that 
#' SNPs are causal. Information on parent case-control status must be included. 
#' 
#' @param child A list object with an FBM.code256 and accompanying FAM and MAP. Must contain case-control status
#' parents in FAM. 
#' @param include Vector of rows to use in regression. Used with cross-validation. Default uses all rows.
#' @return A data.frame with slopes of each regression, standard errors of each slope, t-scores associated with each slope and P-values of each slope.
#' @export 
#' 

GWAX <- function(child, include = bigparallelr::rows_along(child$genotypes)) {
  
  p1_Status <- child$FAM$p1_Status
  p2_Status <- child$FAM$p2_Status
  child_status <- child$FAM$Status
  FBM <- child$genotypes
  
  #Creates a vector of the proxy statuses for the child
  x <- (child_status == 1 | p1_Status == 1 | p2_Status == 1) + 0
  
  GWAS(child, x, include = include)
}

