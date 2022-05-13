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


OpenRds <- function(path){
  rds_file <- paste(path, ".rds", sep = "")
  if (file.exists(rds_file)) {
    snp_attach(rds_file)  
  } else {
    stop("No file with that name.")
  }
} 


