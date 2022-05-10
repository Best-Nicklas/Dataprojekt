#' Description - Helper function to open an Rds file
#' 
#' @param path
#' @return 
#' @example 
#' 

OpenRds <- function(path){
  rds_file <- paste(path, ".rds", sep = "")
  if (file.exists(rds_file)) {
    snp_attach(rds_file)  
  } else {
    stop("No file with that name.")
  }
} 