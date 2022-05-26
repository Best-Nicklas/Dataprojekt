#' Create an .rds file
#' 
#' This function is used to create and save a file-backed matrix with related information.
#' 
#' @param nrow Integer specifying number of rows.
#' @param ncol Integer specifying number of columns.
#' @param path Path where to save the files (DO NOT SPECIFIY FILE EXTENSION). 
#' 
#' @return A list containing a file-backed matrix, a 
#' FAM tibble and a MAP tibble.
#' @details When using this function both a .bk and .rds file will be created at the specified path. Besides 
#' the file-backed matrix two tibbles containing FAM (phenotype information) and MAP (SNP information)
#' will be intialized and saved together with the file-backed matrix in the .rds file. 
#' Both the .bk and the .rds file must always be in the same directory.
#' The function openRds can be used to access these files again (see function reference for openRds). 
#' 
#' @export
#'


createRds <-  function(path, nrow, ncol) {
  if(nrow <= 0 || ncol <= 0 ) stop("nrow and ncol must be 1 or higher")
  
  # Creates a file-backed matrix (.bk file) of the given size with all zeroes.
  G = bigstatsr::FBM.code256(nrow = nrow, # number of rows
                  ncol = ncol, # number of columns
                  code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                  backingfile = path) 
  # creates list containing FAM/MAP info and a pointer to a file_backed matrix 
  obj.bigsnp = list(
    genotypes = G,
    MAP = tibble::tibble(SNP_ID = 1:ncol),
    FAM = tibble::tibble(ID = 1:nrow))
  
  # Saves obj.bigsnp in .rds file. 
  bigsnpr::snp_save(obj.bigsnp)
  
  
  return(obj.bigsnp)
}

