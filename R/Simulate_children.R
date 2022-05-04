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
  
  FBM <- CreateFBM(rows, cols, path, type = NA_real_)
  
  blocks <-  round(seq(0, rows, length = n_blocks + 1))
  for (i in 1:(length(blocks)-1)) {  
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    FBM[b_start:b_end, ] <- truer_gen(p1$genotypes[b_start:b_end, ], 
                                      p2$genotypes[b_start:b_end, ], 
                                      seed)
  }
  
  MAPFAM <- calculate_MAPFAM(FBM, disease, n_blocks, seed)
  
  family_info <- tibble(p1_status = p1$FAM$Status, p2_status = p2$FAM$Status)
  
  children <- list(genotypes = FBM, FAM = MAPFAM$FAM, MAP = MAPFAM$MAP, family_info = family_info)
  snp_save(children)
  
  return(children)
}