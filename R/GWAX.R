#' Description - GWAX
#' 
#' @param child
#' @return 
#' @example 
#' 

GWAX <- function(child) {
  
  n <- length(child$MAP$SNP_ID)
  p1_status <- child$family_info$p1_status
  p2_status <- child$family_info$p2_status
  child_status <- child$FAM$Status
  FBM <- child$genotypes
  
  #Creates a vector of the proxy statuses for the child
  x <- ifelse(child_status == 1 | p1_status == 1 | p2_status == 1, 1, 0)
  
  GWAS(child, x)
}
