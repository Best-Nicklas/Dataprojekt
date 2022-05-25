#' Open Rds
#' 
#' This function is used to open a saved .rds file. If the path do not exists 
#' the function returns a error message. 
#' 
#' @param path Path to file (DO NOT INCLUDE FILE EXTENSION).
#' @return returns the list saved in .rds file. 
#' @export
#' 

OpenRds <- function(path){
  rds_file <- paste(path, ".rds", sep = "")
  if (file.exists(rds_file)) {
    bigsnpr::snp_attach(rds_file)  
  } else {
    stop("No file with that name.")
  }
} 