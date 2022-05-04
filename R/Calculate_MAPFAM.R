#' Description - Function to simulate FAM and MAP for a dataframe with persons´
#' genotypes
#' 
#' @param FBM
#' @param disease
#' @param n_blocks
#' @param seed
#' @return 
#' @example 
#' 


calculate_MAPFAM <- function(FBM, disease, n_blocks, seed) {
  if (!is.null(seed)) {set.seed(seed)}
  
  rows <- nrow(FBM)
  cols <- ncol(FBM)
  
  MAF <- disease$MAF
  beta <- disease$BETA
  causal <- disease$CAUSAL
  h2 <- disease$H2
  prevalens <- disease$PREVALENS
  
  mu <- 2 * MAF * causal
  sigma <- sqrt((2 * MAF * causal)*(1 - MAF * causal))
  sigma[sigma == 0] <- 1
  
  liab_g <- numeric(rows)
  
  blocks <-  round(seq(0, rows, length = n_blocks + 1))
  for (i in 1:(length(blocks) - 1)) {  
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    liab_g[b_start:b_end] <- sweep(sweep(FBM[b_start:b_end, ], 
                                         MARGIN = 2, 
                                         STATS = mu, 
                                         FUN = "-"), 
                                   MARGIN = 2, 
                                   STATS = sigma, 
                                   FUN = "/") %*% beta 
  }
  
  liab_e <- rnorm(rows, 0, sqrt(1 - h2))
  liab_full <- liab_e + liab_g
  threshold <- qnorm(1 - prevalens)
  
  FAM <- tibble(ID = 1:rows, 
                Full_Liability = liab_full, 
                Genetic_Liability = liab_g, 
                Status = ifelse(liab_full > threshold, 1, 0))
  MAP <- tibble(SNP_ID = 1:cols, 
                MAF = disease$MAF, 
                BETA = disease$BETA, 
                CAUSAL = disease$CAUSAL,
                H2 = disease$H2)
  
  return(list(FAM = FAM, MAP = MAP))
}