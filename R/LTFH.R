#' Perform GWAS on posterior mean genetic liabilities
#' 
#' This function finds all the unique status configurations in the data and 
#' calculates the mean posterior genetic liability for each of these unique 
#' configurations using gibbs sampling. The calculated posterior mean genetic 
#' liabilities are then matched to each person in the dataset and used to perform GWAS.
#' 
#' @param child A list object with an FBM.code256 and accompanying FAM and MAP.
#' @param prevalence The likelihood of having the disease in the population
#' @param h2 Heritability parameter.
#' @return The function returns a list, where the first entrance is the GWAS 
#' output and the second entrance are the calculated posterior mean genetic liabilities.
#' @export
#' 

LTFH <- function(child, prevalence, h2) {
  child_configs <- purrr::pmap(list(child$FAM$Status,child$FAM$p1_Status, child$FAM$p2_Status, child$FAM$sibs_Status), 
                        function(x1, x2, x3, x4) 
                          paste(toString(x1),
                                toString(max(x2, x3)),
                                toString(min(x2, x3)),
                                strrep(1, sum(x4)),
                                strrep(0,length(x4) - sum(x4)),
                                sep = "")) %>% do.call("rbind", .)
  
  #Finds all the unique status configurations that the children have
  unique_configs <- unique(child_configs)
  
  #Calculates the the mean posterior genetic liability for each unique config
  config_liabs <- vector()
  for (config in unique_configs) {
    gibb_input <- as.numeric(strsplit(config,"")[[1]])
    config_liabs[config] <- gibbs_sampler(gibb_input, 100, covmatrix(h2 = h2, n_sib = length(gibb_input) - 3), prevalence)
  }
  
  n  <- nrow(child$genotypes)
  gen_liabs <- numeric(n)
  
  #Matches config of child to calculated mean posterior genetic liability
  for (i in 1:n) {
    gen_liabs[i] <- config_liabs[child_configs[i]]
  }
  
  return(list(GWAS_Data = GWAS(child, gen_liabs), Posterior_Mean_Genetic_Liability = gen_liabs))
}