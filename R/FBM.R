#' Description - Function to create File-back matrix
#' 
#' @param nrow 
#' @param ncol 
#' @param path 
#' @param type 
#' @return An empty filebacked matrix 
#' @example 
#'
#'


# FBM file
CreateFBM <- function(nrow, ncol, path, type = NA_integer_){
  return(FBM.code256(nrow = nrow, # number of rows
                     ncol = ncol, # number of columns
                     code = c(0L, 1L, 2L, rep(type, 256 - 3)),
                     backingfile = path)) 
  
}


# Denne funktion tjener ikke rigtigt noget formål. I omdøber bare FBM.code256, hvilket jeg ikke kan se nogen grund til.
# Hvis i gerne vil have en funktion, som i kan indlæse eller lave bigsnp filer med, så synes jeg, at i skal lave en funktion 
# i stil med det vi allerede har snakket om for noget tid siden.

#i pseudo kode halløj
get_FBM = function(path, ...) { #og diverse andre input
  
  if (file.exists(path)) {
    snp_attach(path)
  } else {
    G = FBM.code256(...)
    obj.bigsnp = list(
      genotypes = G,
      map = tibble(...),
      fam = tibble(...)
    )
    
    snp_save(obj.bigsnp)
    return(obj.bigsnp)
  }
}

# dette burde fikse de problemer, som i har snakket om, hvor i ikke har adgang til jeres FBMs, hvis der sker en fejl.
