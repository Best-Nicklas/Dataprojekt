#' Description - 
#' 
#' @param child
#' @return 
#' @example 
#' 

LTFH <- function(child) {
  h2 <- child$MAP$H2
  #Finds all the unique status configurations that the children have
  unique_configs <- unique(child$family_info$config)
  
  #Calculates the the mean posterior genetic liability for each unique config
  config_liabs <- vector()
  for (config in unique_configs) {
    gibb_input <- as.numeric(strsplit(config,"")[[1]])
    config_liabs[config] <- gibbs_sampler(gibb_input, 20000, 1000, covmatrix(h2 = 0.8, sib = length(gibb_input) - 3))
  }
  
  #Loads all configs for children 
  child_configs <- child$family_info$config
  n  <- nrow(child$genotypes)
  gen_liabs <- numeric(n)
  
  #Matches config of child to calculated mean posterior genetic liability
  for (i in 1:n) {
    gen_liabs[i] <- config_liabs[child_configs[i]]
  }
  
  return(gen_liabs)
}
