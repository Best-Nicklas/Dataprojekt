#' Description - Simulation child 
#' 
#' @param p1
#' @param p2
#' @param disease
#' @param path
#' @param n_blocks
#' @param seed
#' @return 
#' @example 
#' 


Simulate_children <- function(p1, p2, disease, path, n_blocks = 20, overwrite = T, seed = NULL) {
  rows <- nrow(p1$genotypes)
  cols <- ncol(p1$genotypes)
  
  #Initializes an empty FBM or opens an existing one 
  FBM <- verifyRds(path, overwrite, rows, cols)
  
  #Calculates blocks from specified number of blocks n_blocks
  blocks <-  round(seq(0, rows, length = n_blocks + 1))
  
  #Creates children in block sizes 
  for (i in 1:(length(blocks)-1)) {  
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    FBM$genotypes[b_start:b_end, ] <- truer_gen(p1$genotypes[b_start:b_end, ], 
                                                p2$genotypes[b_start:b_end, ], 
                                                seed)
  }
  
  #Calculates MAF and FAM info for the children
  MAPFAM <- calculate_MAPFAM(FBM$genotypes, disease, n_blocks, seed)
  FBM$MAP <- MAPFAM$MAP
  FBM$FAM <- MAPFAM$FAM
  
  #Creates family info for the children
  FBM$FAM$p1_status <- p1$FAM$Status
  FBM$FAM$p2_status <- p2$FAM$Status
  
  return(FBM)
}