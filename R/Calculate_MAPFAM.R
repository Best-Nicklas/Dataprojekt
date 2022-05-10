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


calculate_MAPFAM <- function(FBM, disease, n_blocks = 20, seed) {
  if (!is.null(seed)) {set.seed(seed)}
  
  #Initializes parameters
  rows <- nrow(FBM)
  cols <- ncol(FBM)
  
  MAF <- disease$MAF
  beta <- disease$BETA
  causal <- disease$CAUSAL
  h2 <- disease$H2
  prevalens <- disease$PREVALENS
  
  #calculates mu and sigma used for normalization
  mu <- 2 * MAF * causal
  sigma <- sqrt((2 * MAF * causal)*(1 - MAF * causal))
  sigma[sigma == 0] <- 1
  
  #Determines the blocks based on specificed number n_blocks
  blocks <-  round(seq(0, rows, length = n_blocks + 1))
  
  #Calculates genetic liablities in block sizes 
  liab_g <- future_lapply(1:(length(blocks) - 1), function(i) {
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    sweep(sweep(FBM[b_start:b_end, ], 
                MARGIN = 2, 
                STATS = mu, 
                FUN = "-"), 
          MARGIN = 2, 
          STATS = sigma, 
          FUN = "/") %*% beta
    
  }, future.seed = T) %>% do.call("rbind", .) %>% as.numeric()
  
  threshold <- qnorm(1 - prevalens)
  
  #Creates tibbles that store info disease and liability
  FAM <- tibble(ID = 1:rows, 
                Genetic_Liability = liab_g,
                Full_Liability = rnorm(rows, 0, sqrt(1 - h2)) + Genetic_Liability, 
                Status = (Full_Liability > threshold) + 0)
  
  MAP <- tibble(SNP_ID = 1:cols, 
                MAF = disease$MAF, 
                BETA = disease$BETA, 
                CAUSAL = disease$CAUSAL,
                H2 = disease$H2)
  
  return(list(FAM = FAM, MAP = MAP))
}