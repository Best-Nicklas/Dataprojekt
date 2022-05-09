#' Description - Helper function to calculate siblings liability 
#' 
#' @param n_sibs
#' @param p1
#' @param p2
#' @param disease
#' @return 
#' @example 
#' 
# hvis i gerne vil have n_sib som argument, så burde i også bruge den lidt mere.
# Nu hedder funktionen jo sibs_liab_calc. så hvorfor ikke lave en funktion, som laver n_sibs søskende og returnere deres relevante information?
# på den måde kan i bare kalde den efter i har lavet jeres første barn eller bare returnere genotyperne for den sidst udregnet "søskende" og lave
# den til jeres genotyped barn.
sibs_liab_calc <- function(n_sibs, p1, p2, disease) { 
  if (n_sibs == 0) {
    return(NULL)
  }
  
  MAF <- disease$MAF
  beta <- disease$BETA
  causal <- disease$CAUSAL
  
  #Generates children from p1 and p2
  temp <- truer_gen(p1, p2, NULL)
  
  #Calculates mu and sigma used for normalization
  mu <- 2 * MAF * causal
  sigma <- sqrt((2 * MAF * causal)*(1 - MAF * causal))
  sigma[sigma == 0] <- 1
  
  c <- sum(causal)
  
  #Calculates the genetic liabilities for sibs
  # i stedet for at kalde dette udtryk hele tiden, så lav en funktion der hedder normalize_genotypes eller noget, som tager 
  # en mat som input og tilsvarende mu, sigma og beta.
  liab_g <- c(sweep(sweep(temp, 
                          MARGIN = 2, 
                          STATS = mu, 
                          FUN = "-"), 
                    MARGIN = 2, 
                    STATS = sigma, 
                    FUN = "/") %*% beta)
  return(liab_g)  
}

