#' Open an .rds file
#' 
#' This function is used to open a saved .rds file. If the path does not exist 
#' the function returns a error message. 
#' 
#' @param path Path to file (DO NOT INCLUDE FILE EXTENSION).
#' @return returns the rds object saved in .rds file. 
#' @export
#' 

OpenRds <- function(path){
  # Format correct file name
  rds_file <- paste(path, ".rds", sep = "")
  
  # Open file or give error
  if (file.exists(rds_file)) {
    bigsnpr::snp_attach(rds_file)  
  } else stop("No file with that name.")
} 