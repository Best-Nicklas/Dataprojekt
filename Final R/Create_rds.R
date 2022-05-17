#' Description - Create Rds
#' @param path Path to file
#' @param nrow Total rows
#' @param ncol Total columns 
#' @return 
#' @example 
#' 

createRds <-  function(path, nrow, ncol) {
  G = FBM.code256(nrow = nrow, # number of rows
                  ncol = ncol, # number of columns
                  code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                  backingfile = path) 
  
  obj.bigsnp = list(
    genotypes = G,
    MAP = tibble(SNP_ID = 1:ncol),
    FAM = tibble(ID = 1:nrow))
  snp_save(obj.bigsnp)
  
  
  return(obj.bigsnp)
}