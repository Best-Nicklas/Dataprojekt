#' Description - Verify Rds
#' @param path Path to file
#' @param overwrite Boolean value 
#' @param nrow Total rows
#' @param ncol Total columns
#' @return 
#' @example 
#' 

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
