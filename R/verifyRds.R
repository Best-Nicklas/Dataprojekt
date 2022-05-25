#' Verification of .rds file
#' 
#' This is an internal helper function used to verify if the path is valid. If the path is already in use, 
#' it checks if the already exiting .rds file has the same dimensions and overwrites it (if overwrite TRUE). 
#' If the dimensions are wrong it stops and returns a error message. 
#' 
#' @param path Path where to save and find the file. 
#' @param overwrite Boolean value specifying whether to overwrite exisiting files.  
#' @param nrow Integer specifying the amount of rows.
#' @param ncol Integer specifying the amount of columns. 
#' @keywords internal
#' @export


verifyRds <- function(path, overwrite, nrow, ncol) {
  rds_file <- paste(path, ".rds", sep = "")
  if (file.exists(rds_file)) {
    if (overwrite) {
      FBM <- OpenRds(path)
      if (nrow(FBM$genotypes) != nrow | ncol(FBM$genotypes) != ncol) stop("An .rds file exists with same name but wrong dimensions. 
                                                           Cannot overwrite. Please choose another filename or change dimensions.")
    } else {
      stop("An .rds file of this name already exists. 
            Please allow overwrite or choose another filename.")
    }
  } else {
    FBM <- createRds(path, nrow, ncol)
  }
  return(FBM)
} 
