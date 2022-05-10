#' Description - Function to create or load filebacked matrices and related information from or to a rds file
#' 
#' @param nrow 
#' @param ncol 
#' @param path 
#' @param type 
#' @return An rds object containing a filebacked matrix and as well as FAM and MAP information
#' @example 
#'
#'


verifyRds <- function(path, overwrite, nrow, ncol) {
  rds_file <- paste(path, ".rds", sep = "")
  if (file.exists(rds_file)) {
    if (overwrite) {
      FBM <- OpenRds(path)
      if (nrow(FBM$genotypes) != nrow | ncol(FBM$genotypes) != ncol) stop("An .rds file exists with same name but wrong dimensions. 
                                                           Cannot overwrite. Please choose another filename or change dimensions")
    } else {
      stop("An .rds file of this name already exists. 
            Please allow overwrite or choose another filename")
    }
  } else {
    FBM <- createRds(path, nrow, ncol)
  }
  return(FBM)
} 

