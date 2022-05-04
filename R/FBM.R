#' Description - Function to create File-back matrix
#' 
#' @param nrow 
#' @param ncol 
#' @param path 
#' @param type 
#' @return The File-back matrxi
#' @example 
#'
#'


# FBM file
CreateFBM <- function(nrow, ncol, path, type = NA_integer_){
  return(FBM.code256(nrow = nrow, # number of rows
                     ncol = ncol, # number of columns
                     code = c(0L, 1L, 2L, rep(type, 256 - 3)),
                     backingfile = path)) 
  
}
