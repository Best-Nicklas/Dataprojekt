#' Function to simulate FAM and MAP files 
#'
#' @param dataframe 
#' @param MAF 
#' @param causal_SNP
#' @param h2
#' @param seed 
#' @return 
#' @examples   
#' Simulate_FAM_MAP()

Simulate_FAM_MAP <- function(dataframe, MAF, causal_SNP , h2, seed = NULL) {
  if (!is.null(seed)) {set.seed(seed)}
  
  rows <- dim(dataframe)[1]
  cols <- dim(dataframe)[2]
  
  if (cols != length(MAF)) stop("Number of columns in 
                                 dataframe does not match 
                                 length of MAF vector")
  
  c <- length(causal_SNP)
  beta <- numeric(cols)
  beta[causal_SNP] <- rnorm(length(causal_SNP), 0, sqrt(h2/c))
  
  mu <- 2*MAF[causal_SNP]
  sigma <- sqrt((2*MAF[causal_SNP])*(1-MAF[causal_SNP]))
  
  liab_g <- numeric(rows)
  liab_g <- sweep(sweep(dataframe[,causal_SNP], MARGIN = 2, STATS = mu, FUN = "-"), MARGIN = 2, STATS = sigma, FUN = "/") %*%    
    beta[causal_SNP]
  
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
                Causal = ifelse(beta != 0, 1, 0))
  
  MAP_and_FAM <- list(genotypes = dataframe, fam = FAM, map = MAP)
  snp_save(MAP_and_FAM)
}

#' Helper function to open Map and Fam file
#'
#' @param path The path where we have saved the FBM 
#' @return 
#' @examples   
#' OpenRds("../Data/Simulate")

OpenRds <- function(path){
  snp_attach(path)  
}
