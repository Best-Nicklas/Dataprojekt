#' Description - 
#' 
#' @param n
#' @param disease
#' @param filename
#' @param n_blocks
#' @param seed
#' @return 
#' @example 
#' 


simulate_people <- function(n, disease, file_name, n_blocks = 20, seed = NULL) {
  rows <- n
  cols <- disease$N_SNP
  
  # bk_file <- paste(file_name, ".bk", sep = "")
  # rds_file <- paste(file_name, ".rds", sep = "")
  # 
  # if (file.exists(bk_file) & file.exists(rds_file) & overwrite) {
  #   file.remove(bk_file)
  #   file.remove(rds_file)
  # }
  # else if (file.exists(paste(file_name, ".bk", sep = "")) & overwrite) {
  #   file.remove(bk_file)
  # }
  
  #Initializes empty FBM 
  FBM <- CreateFBM(rows, cols, file_name) ### erstat med en opdateret funktion, så i ikke mister jeres FBM.
  
  #Fills empty FBM out with genotypes
  PopulateFBM(FBM, disease$MAF, seed = seed) #eventuelt sig: null_catcher = populateFBM() for at fange nulls, så console ikke bliver spammet
  
  #Calculates MAP and FAM info for the genotypes
  MAPFAM <- calculate_MAPFAM(FBM, disease, n_blocks, seed)
  
  #Saves info in Rds file
  all_people <- list(genotypes = FBM, FAM = MAPFAM$FAM, MAP = MAPFAM$MAP)
  snp_save(all_people)
  
  return(all_people)
}
