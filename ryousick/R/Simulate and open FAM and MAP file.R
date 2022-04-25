#' Function to simulate FAM and MAP files 
#'
#' @param dataframe 
#' @param MAF 
#' @param causal_SNPs
#' @param h2
#' @param seed 
#' @return 
#' @examples   

Simulate_FAM_MAP <- function(dataframe, MAF, causal_SNPs , h2, seed = NULL) {
  if (!is.null(seed)) {set.seed(seed)}
  
  rows <- dim(dataframe)[1]
  cols <- dim(dataframe)[2]
  
  if (cols != length(MAF)) stop("Number of columns in 
                                 dataframe does not match 
                                 length of MAF vector")
  c <- sum(causal_SNP == 1)
  beta <- ifelse(causal_SNP == 1, rnorm(cols, 0, sqrt(h2/c)), 0)
  
  mu <- 2*MAF*causal_SNP 
  sigma <- sqrt((2*MAF*causal_SNP)*(1-MAF*causal_SNP))
  sigma[sigma == 0] <- 1 
  liab_g <- numeric(rows)
  
  b_size <- cols/20
  for (i in 0:((cols/b_size)-1)) {  
    b_start <- 1 + i*b_size
    b_end <- (1+i)*b_size
    liab_g[b_start:b_end] <- sweep(sweep(dataframe[b_start:b_end, ], MARGIN = 2, STATS = mu, FUN = "-"), MARGIN = 2, STATS = sigma, FUN = "/") %*% beta }
  
  liab_e <- rnorm(dim(dataframe)[1], 0, sqrt(1-h2))
  liab_full <- liab_e + liab_g
  threshold <- qnorm(0.95, 0, 1)
  
  FAM <- tibble(ID = 1:rows, 
                Full_Liability = liab_full, 
                Genetic_Liability = liab_g, 
                Status = ifelse(liab_full > threshold, 1, 0))
  MAP <- tibble(SNP_ID = 1:cols, 
                Minor_Allele_Frequency = MAF, 
                BETA = beta, 
                Causal = causal_SNP)
  
  MAP_and_FAM <- list(genotypes = dataframe, fam = FAM, map = MAP)
  snp_save(MAP_and_FAM)
}

#' Helper function to open Map and Fam file
#'
#' @param path The path where we have saved the FBM 
#' @return 
#' @examples   
#' OpenRds("../Data/Simulation_1")

OpenRds <- function(path){
  snp_attach(path)  
}
