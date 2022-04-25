#' Create a file back matrix
#'
#' @param nrow How many persons we want to simulate
#' @param ncol How many SNP´s each persons has
#' @param path The path were we want to save the FBM
#' @return A file back matrix with lenght nrow and ncol with a specific path
#' @examples   
#' CreateFBM(1000, 1000, "../Data/Simulate")
#' CreateFBM(10000, 10000, "../Data/Simulate2")

CreateFBM <- function(nrow, ncol, path, type = NA_integer_){
  return(FBM.code256(nrow = nrow, # number of rows
                     ncol = ncol, # number of columns
                     code = c(0L, 1L, 2L, rep(type, 256 - 3)),
                     backingfile = path)) 
  
}
