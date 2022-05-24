#' Description - 
#' 
#' Open Rds
#' @param path Path to file
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