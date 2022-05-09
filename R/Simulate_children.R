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
  cols <- ncol(p1$genotypes) #overvej at simulér forældre i blocks, e.g. 1000 ad gangen, så kan i undgå at lave to ekstra fbms
  
  #Initializes an empty FBM
  FBM <- CreateFBM(rows, cols, path, type = NA_real_) # erstat med noget der gemmer FBM med dummy map og fam, eller læser fra en path
  
  #Calculates blocks from specified number of blocks n_blocks
  blocks <-  round(seq(0, rows, length = n_blocks + 1))

  #Creates children in block sizes 
  for (i in 1:(length(blocks)-1)) {  # lav til future.apply i stedet. brug seq_along
    
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    FBM[b_start:b_end, ] <- truer_gen(p1$genotypes[b_start:b_end, ], 
                                      p2$genotypes[b_start:b_end, ], 
                                      seed)
    NULL #bare returnér null, ellers vil den vist returnere det sidste i har udregnet, i.e. jeres genotyper. 
    # undgå generelt at flytte data i mellem parent og child processses når i parallelisere. 
    # det kan tage så lang tid, at i mister alle fordele ved at parallelisere, fordi i bare flytter data rundt.
    # I stedet kan i indlæse data i jeres paralleliseret loop og skrive dens resultater til disk.
    # det er ofte hurtigere end at sende mellem parallele processer.
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