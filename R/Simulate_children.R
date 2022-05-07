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


Simulate_children <- function(p1, p2, disease, path, n_blocks = 20, seed = NULL) {
  rows <- nrow(p1$genotypes)
  cols <- ncol(p1$genotypes)
  
  #Initializes an empty FBM
  FBM <- CreateFBM(rows, cols, path, type = NA_real_)
  
  #Calculates blocks from specified number of blocks n_blocks
  blocks <-  round(seq(0, rows, length = n_blocks + 1))

  #Creates children in block sizes 
  for (i in 1:(length(blocks)-1)) {  
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    FBM[b_start:b_end, ] <- truer_gen(p1$genotypes[b_start:b_end, ], 
                                      p2$genotypes[b_start:b_end, ], 
                                      seed)
  }
  
  #Calculates MAF and FAM info for the children
  MAPFAM <- calculate_MAPFAM(FBM, disease, n_blocks, seed)
  
  #Creates family info for the children
  family_info <- tibble(p1_status = p1$FAM$Status, p2_status = p2$FAM$Status)
  
  #Saves all info as Rds file
  children <- list(genotypes = FBM, FAM = MAPFAM$FAM, MAP = MAPFAM$MAP, family_info = family_info)
  snp_save(children)
  
  return(children)
}