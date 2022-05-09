#' Description - 
#' 
#' @param child
#' @return 
#' @example 
#' 

# jeg vil mene, at i skal gentænke jeres input til denne funktion.
# forestil jer, at i har en phenotype fil, i.e. i har id på alle (genotyped) personer, deres status, forældrenes status, og eventuelle søskende status
# Hvordan kan i komme fra sådan en fil ned til noget i kan bruge jeres gibbs sampler på ?
# I skal med andre ord kunne udtrække alt den relevante information fra den fil om konfigurationer osv.
# gerne prøv med en løsning der bygger på tibbles, dplyr og lignende.
# start fx med at konstruere jeres configs fra jeres phenotype fil, og lad siblings (som kan variere i antal) være en liste indgang.

LTFH <- function(child) { # hav h2 som arguement
  h2 <- child$MAP$H2
  #Finds all the unique status configurations that the children have
  unique_configs <- unique(child$family_info$config)
  
  #Calculates the the mean posterior genetic liability for each unique config
  config_liabs <- vector()
  for (config in unique_configs) {
    gibb_input <- as.numeric(strsplit(config,"")[[1]])
    config_liabs[config] <- gibbs_sampler(gibb_input, 20000, 1000, covmatrix(h2 = 0.8, sib = length(gibb_input) - 3)) # i skal ikke hardcode h2
  } #hvordan afgører i om i har et "præcist" nok estimat? sem
  
  #Loads all configs for children 
  child_configs <- child$family_info$config
  n  <- nrow(child$genotypes)
  gen_liabs <- numeric(n)
  
  #Matches config of child to calculated mean posterior genetic liability
  for (i in 1:n) {
    gen_liabs[i] <- config_liabs[child_configs[i]]
  }
  
  # Hvis i bruger tibbles, så kan i have en søjle med configs i jeres phenotype fil og i jeres estimat tibble, og så matche på configs.
  # det burde være meget hurtigt og let, og så undgår i unødvendige loops
  
  return(gen_liabs) # returnér gerne person ids med deres estimat. ellers bliver det (potentielt) bøvlet at holde styr på folk.
}
