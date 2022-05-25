#' Simulation of genotypes with no family
#' 
#' This function is used to simulate genotypes individuals with no family and accompanying phenotype data. 
#' 
#' @param n Integer specifying amount of gentypes/indivduals to simulate. 
#' @param disease A list with all the disease parameters.
#' @param path Path to where .rds file should be saved, or where one is stored if overwriting existing .rds file (DO NOT SPECIFY FILE EXTENSION).
#' @param overwrite Boolean value used to determine if existing .rds file with specified name should be overwritten (Default value TRUE).
#' @param n_blocks Integer used to determine number of blocks to run simulation in (Default value is 20). Set higher if running into memory issues.
#' @details Simulation can be performed using paralization if a paralization plan has been set prior to excution in the global environment. 
#' @return Returns list object containing an FMB.code256 with genotypes, MAF object containing information on SNPs and
#' FAM object containing phenotype information on genotypes
#' 

sim_genotypes_no_family <- function(n, disease, path, overwrite = T, n_blocks = 20) {
  
  # Load disease information
  cols <- disease$N_SNP
  MAF <- disease$MAF
  beta <- disease$BETA
  causal <- disease$CAUSAL
  prevalence <- disease$PREVALENCE
  h2 <- disease$H2
  
  # Calculate normalization constants 
  norm_const <- calc_normalization_consts(MAF, causal)
  mu <- norm_const$mu
  sigma <- norm_const$sigma
  
  #Checks if a FBM with the given name and dimensions exists, else creates one
  FBM <- verifyRds(path, overwrite, n, cols)
  
  # Prepare correct block indexes
  blocks <-  round(seq(0, n, length = n_blocks + 1))
  
  # Inserts values into FBM and calculate genetic liabilities in block sizes
  g_liabs <- future.apply::future_lapply(1:(length(blocks) - 1), function(i) {  
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    FBM_temp <- matrix(rbinom(cols * b_size, 2, MAF),
                       nrow = b_size,
                       byrow = T)
    
    FBM$genotypes[b_start:b_end, ] <- FBM_temp
    
    g_block <- calc_gliab(FBM_temp, beta, mu, sigma)
    
  }, future.seed = T) %>% do.call("rbind", .) %>% as.numeric()
  
  # Saves liability and status information as well as SNP information in Rds
  threshold <- qnorm(prevalence, lower.tail = F)
  FBM$FAM$Genetic_Liability <- g_liabs
  FBM$FAM$Full_Liability <- g_liabs + rnorm(n, 0, sqrt(1 - h2))
  FBM$FAM$Status <- (FBM$FAM$Full_Liability > threshold) + 0
  
  FBM$MAP$MAF <- MAF
  FBM$MAP$BETA <- beta
  
  return(FBM)
}
