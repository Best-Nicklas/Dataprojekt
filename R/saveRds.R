#' Saving the work done to an rds object
#' 
#' This function is used to save the changes made to an rds object during an R session.
#' 
#' @param rds.obj The rds object to save.
#' @export


saveRds <- function(rds.obj) {
  # Wrapper to conform with RyouSick name conventions
  bigsnpr::snp_save(rds.obj)
} 
