#' Description
#'
#' @param MAF 
#' @param dataframe
#' @param seed 
#' @return 
#' @examples 

PopulateFBM <- function(MAF, dataframe, seed = NULL) {
  if (!is.null(seed)) {set.seed(seed)}
  if (dim(dataframe)[2] != length(MAF)) stop("Number of columns in 
                                              dataframe does not match 
                                              length of MAF vector")
  n <- dim(dataframe)[1]
  b_size <- n/20     
  
  for (i in 0:((n/b_size)-1)) {  
    b_start <- 1 + i*b_size
    b_end <- (1+i)*b_size
    dataframe[b_start:b_end,] <- matrix(nrow=b_size,
                                        rbinom(n*b_size,2,MAF), 
                                        byrow=TRUE)
  }
}
