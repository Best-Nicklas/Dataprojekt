#' Description - Function to populate the FBM
#' 
#' @param MAF
#' @param dataframe
#' @param seed
#' @example 
#' 

PopulateFBM <- function(FBM, MAF, n_blocks = 20, seed = NULL) {
  #i kan spare jer selv for meget bøvl ved at lave denne slags tjek. stopper kode tidligt, hvis simple ting ikke er opfyldt.
  if (ncol(FBM) != length(MAF)) stop("Number of columns in 
                                              dataframe does not match 
                                              length of MAF vector")
  
  
  if (!is.null(seed)) {set.seed(seed)}
  
  # Finds the correct blocks using the n_blocks parameter
  blocks <-  round(seq(0, nrow(FBM), length = n_blocks + 1))
  
  # Inserts values into FBM in block sizes
  for (i in 1:(length(blocks) - 1)) {   # i kan med fordel lave dette loop til et future.lapply loop med null return. brug seq_along
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    FBM[b_start:b_end, ] <- matrix(nrow = b_size,
                                   rbinom(ncol(FBM) * b_size, 2, MAF), 
                                   byrow = TRUE)
    NULL #returnere bare null nu
  }
  NULL #returnere bare null her. i skal ikke bruge noget return.
}


