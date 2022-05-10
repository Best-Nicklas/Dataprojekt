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


createRds <-  function(path, nrow, ncol) {
  G = FBM.code256(nrow = nrow, # number of rows
                  ncol = ncol, # number of columns
                  code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                  backingfile = path) 
  
  obj.bigsnp = list(
    genotypes = G,
    MAP = tibble(ID = 1:nrow),
    FAM = tibble(SNP_ID = 1:ncol))
  snp_save(obj.bigsnp)
  
  
  return(obj.bigsnp)
}


