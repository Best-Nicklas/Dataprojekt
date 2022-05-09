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
  
  
  #det er bedre kode at lave et loop med fx seq_along(blocks), da i ikke kan løbe ind i underlige problemer som fx 1:0 eller andet besynderligt.
  # seq_along på tomme objekter returnere bare null, i.e. intet loop bliver lavet.
  liab_g <- future_lapply(1:(length(blocks) - 1), function(i) { 
    
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    # I skal passe meget på med at returnere data fra jeres paralleliering.
    # I vil nok helst bare lave udregningerne, og så gemme dem direkte i jeres genotyper
    # I kan godt returnere mindre data, som fx liabilities eller summary info, men ikke hele genotype matricer.
    
    sweep(sweep(FBM[b_start:b_end, ], 
                MARGIN = 2, 
                STATS = mu, 
                FUN = "-"), 
          MARGIN = 2, 
          STATS = sigma, 
          FUN = "/") %*% beta
    
  }, future.seed = T) %>% do.call("rbind", .) %>% as.numeric()
  
  #Adds environmental factor
#  liab_full <- liab_e + liab_g ikke nødvendigt at lave dette objekt. alle udregningerne med denne sker i jeres tibble.
#  liab_e <- rnorm(rows, 0, sqrt(1 - h2)) samme her
  threshold <- qnorm(1 - prevalens)# hvis værdierne er små, så kan dette udtryk godt fejle.
  #overvej i stedet at bruge qnorm(0.05, lower.tail = F). Samme svar værdi, men mere stabil overfor små værdier.
  
  #Creates tibbles that store info disease and liability
  FAM <- tibble(ID = 1:rows, 
                Genetic_Liability = liab_g, 
                Full_Liability = rnorm(rows, 0, sqrt(1 - h2)) + Genetic_Liability, #simuler bare environment her, så gemmer i ikke noget unødvendigt
                Status = (liab_full > threshold) + 0)  #ifelse(liab_full > threshold, 1, 0)) ifelse er kun nødvendig, hvis man vil have andet end 0 og 1
  MAP <- tibble(SNP_ID = 1:cols, 
                MAF = disease$MAF, 
                BETA = disease$BETA, 
                CAUSAL = disease$CAUSAL, # denne info er overflødig, da vi allerede ved hvilke der er kausale udfra betaerne.
                H2 = disease$H2) #her gentager i bare h2 en masse gange. overvej om i ikke kan gemme denne info smartere. fx i navnene på filerne.
  
  return(list(FAM = FAM, MAP = MAP))
}