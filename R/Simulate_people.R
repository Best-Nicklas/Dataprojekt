#' Description - 
#' 
#' @param n
#' @param disease
#' @param filename
#' @param n_blocks
#' @param seed
#' @return 
#' @example 
#' 


simulate_people <- function(n, disease, path, n_blocks = 20, overwrite = T, seed = NULL) {
  rows <- n
  cols <- disease$N_SNP
  
  #Checks if a FBM with the given name and dimensions exists, else creates one
  FBM <- verifyRds(path, overwrite, rows, cols) 
  
  #Fills empty FBM out with genotypes
  PopulateFBM(FBM$genotypes, disease$MAF, seed = seed)
  
  #Calculates MAP and FAM info for the genotypes
  MAPFAM <- calculate_MAPFAM(FBM$genotypes, disease, n_blocks, seed)
  
  #Saves info in Rds file
  FBM$FAM <- MAPFAM$FAM 
  FBM$MAP <- MAPFAM$MAP
  
  return(FBM)
}