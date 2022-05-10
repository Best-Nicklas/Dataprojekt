#' Description - 
#' 
#' @param child
#' @return 
#' @example 
#' 

LTFH <- function(child) {
  h2 <- child$MAP$H2 
  child_configs <- pmap(list(child$FAM$Status,child$FAM$p1_status, child$FAM$p2_status, child$FAM$sibs_status), 
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
    config_liabs[config] <- gibbs_sampler(gibb_input, 100, covmatrix(h2 = 0.8, n_sib = length(gibb_input) - 3))
  }
  
  n  <- nrow(child$genotypes)
  gen_liabs <- numeric(n)
  
  #Matches config of child to calculated mean posterior genetic liability
  for (i in 1:n) {
    gen_liabs[i] <- config_liabs[child_configs[i]]
  }
  
  return(GWAS(child, gen_liabs))
}
