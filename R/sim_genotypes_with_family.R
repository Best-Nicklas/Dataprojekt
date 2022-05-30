#' Simulation of family (Fixed and non fixed)
#' 
#' This function is used to simulate genotypes gentypes for individuals with accompanying phenotype information and 
#' information on the case-control status of parents annd possible sibling.
#' 
#' @param n Integer specifying the number of individuals/genotypes to simulate. 
#' @param disease A list with all the disease parameters.
#' @param path Path to where .rds file should be saved, or where one is stored if overwriting existing .rds file (DO NOT SPECIFY FILE EXTENSION).
#' @param n_sibs Integer value for how many sibling to produce for each genotype or vector containing values for
#' how many sibs to sample from (ex. c(1,4,6) will produce genotypes randomly with 1, 4 or 6 siblings).  
#' @param overwrite Boolean value used to determine if a helper function is allowed to overwrite (Default value is TRUE)
#' @param n_blocks Integer used to determine number of blocks to run simulation in (Default value is 500). Set higher if running into memory 
#' issues such as freezing or crashing. Setting n_blocks higher reduces the memory size of each block, but slightly slows the calculation time.
#' @return Returns list object containing an FMB.code256 with genotypes, MAF object containing information on SNPs and
#' FAM object containing phenotype information on genotypes as well as case-control status of parents and possible siblings.
#' @details Simulating a 100.000x100.000 dataset will take up around 9.76 GB of space. Since the running time depends on a number of variables,
#' such as the parallelization settings, core speed, core amount and sibling configuration, we cannot accurately give an estimation how long the simulation will take.
#' Instead we simply warn the user that simulations might take upwards of multiple hours for large datasets such as a 100.000x100.000. The default n_blocks parameter has been set to 500 as this
#' is the number at which a 100.000x100.000 with 2 siblings for each genotype will use a maximum of 2 GB of RAM for calculating a single block. 
#' Simulation can be performed using parallelization if a parallelization plan has been set prior to execution in the global environment. 
#' WARNING: using parallelization will, with a n_blocks of 500, use up to a maximum of 2 GB of RAM for EACH process 
#' when running a simulation of 100.000x100.000 with 2 siblings for each genotype.
#' @export
#' 

sim_genotypes_with_family <- function(n, disease, path, n_sibs = NULL, overwrite = T, n_blocks = min(n, 500)) {
  if (n <= 0) stop("n must positive")
  if (any(n_sibs < 0)) stop("all sibling values must be positive")
  
  # Load disease information
  cols <- disease$N_SNP
  MAF <- disease$MAF
  beta <- disease$BETA
  causal <- disease$CAUSAL
  prevalence <- disease$PREVALENCE
  h2 <- disease$H2
  
  # Create or find FBM file to fill with child genotypes
  FBM <- verifyRds(path, overwrite, n, cols)
  
  # Create vector with number of sibs for each child
  if (!is.null(n_sibs)) {
    sibs_pr_child <- if (length(n_sibs) == 1) {rep(n_sibs,n)} else 
      sample(n_sibs, n, replace = T)
  }
  
  # Calculate normalization constants 
  norm_const <- calc_normalization_consts(MAF, causal)
  mu <- norm_const$mu
  sigma <- norm_const$sigma
  
  #prepare for block calculations of genetic liabilities 
  blocks <-  round(seq(0, n, length = n_blocks + 1))
  g_liabs <- future.apply::future_lapply(1:(length(blocks) - 1), function(i) {  
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    # simulate parent 1 genotypes and genetic liabilities 
    p1 <- matrix(rbinom(cols * b_size, 2, MAF),
                 nrow = b_size,
                 byrow = T)
    p1_gliab <- calc_gliab(p1, beta, mu, sigma)
    
    # simulate parent 2 genotypes and genetic liabilities 
    p2 <- matrix(rbinom(cols * b_size, 2, MAF),
                 nrow = b_size,
                 byrow = T)
    p2_gliab <- calc_gliab(p2, beta, mu, sigma)
    
    # Simulate child genotypes and genetic liabilities - store genotypes in FBM
    child <- child_gen(p1, p2)
    FBM$genotypes[b_start:b_end, ] <- child
    child_gliab <- calc_gliab(child, beta, mu, sigma)
    
    # Generate sibs for each child and calculate their genetic liabilities
    sibs_gliab <- vector(mode = "list", b_size)
    if (!is.null(n_sibs)) {
      for (s in 1:b_size) {
        sibs_pos <- b_start + s - 1
        if (sibs_pr_child[sibs_pos] != 0) {
          sibs <- child_gen(matrix(rep(p1[s, ], sibs_pr_child[sibs_pos]), 
                                   sibs_pr_child[sibs_pos], cols, byrow = T),
                            matrix(rep(p2[s, ], sibs_pr_child[sibs_pos]), 
                                   sibs_pr_child[sibs_pos], cols, byrow = T))
          sibs_gliab[[s]] <- calc_gliab(sibs, beta, mu, sigma)  
        } else sibs_gliab[[s]] <- NULL
        
      }
    }
    
    tibble::tibble(child_gliab, p1_gliab, p2_gliab, sibs_gliab)
    
  }, future.seed = T) %>% do.call(dplyr::bind_rows, .)
  
  # Calculate full liabilities/status and insert in rds file object
  threshold <- qnorm(prevalence, lower.tail = F)
  FBM$FAM$Genetic_Liability <- g_liabs$child_gliab
  FBM$FAM$Full_Liability <- g_liabs$child_gliab + rnorm(n, 0, sqrt(1 - h2))
  FBM$FAM$Status <- (FBM$FAM$Full_Liability > threshold) + 0
  FBM$FAM$p1_Status <- (g_liabs$p1_gliab + rnorm(n, 0, sqrt(1 - h2)) > threshold) + 0
  FBM$FAM$p2_Status <- (g_liabs$p2_gliab + rnorm(n, 0, sqrt(1 - h2)) > threshold) + 0
  
  FBM$MAP$MAF <- MAF
  FBM$MAP$BETA  <- beta
  
  if (!is.null(n_sibs)) {
    sibs_Full_Liability <- purrr::map(g_liabs$sibs_gliab, .f = ~ 
                                        {if ( is.null(.x)) NULL 
                                          else .x + rnorm(length(.x), 0, sqrt(1 - h2))})
    FBM$FAM$sibs_Status <- purrr::map(sibs_Full_Liability, .f = ~ 
                                        {if ( is.null(.x)) NULL
                                          else(.x > threshold) + 0})
  }
  bigsnpr::snp_save(FBM)
  return(FBM)
  
  
}