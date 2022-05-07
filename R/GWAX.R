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
  x <- numeric(n)
  for (i in 1:n) {
    if (child_status[i] == 1 | p1_status[i] == 1 | p2_status[i] == 1) {
      x[i] <- 1  
    }
    else {
      x[i] <- 0
    }
  }
  
  #Uses function from bigSNPr package perform regression on the FBM
  regr <- big_univLinReg(FBM, x)
  regr$p.value <- predict(regr, log10 = FALSE)
  
  return(data.frame(regr)) 
}
