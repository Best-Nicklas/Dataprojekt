library(bigsnpr)

G2 = FBM.code256(nrow = 1e5, # number of rows
                ncol = 1e5, # number of columns
                code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)), # coding of FileBackedMatrix, here 0,1,2 and NA's.
                backingfile = "test_fbm") # where to save the matrix on disk?

#coding indicates which values can be stored in the FBM. specified with a vector of length 256. 

#access genotypes as usual:
G2[1:10,1:10]

# NEVER (!!!!) use G[,]. It attempts to load ALL of the file backed matrix into memory.
# most efficient approach is to fill out columns. This is due to how objects are accessed. 
# But filling out rows first might be fine too. 

# adhoc way to fill out the matrix
for(i in 1:10){
  #fill out G in a loop with some block size
}

#filling with some dummy vars.
G2[,1:10] <- matrix(0:2, ncol = 10, nrow = 1e5)
G2[1:10,1:10]

#some functions in this package can very efficiently perform calculations on FBMs, eg:
tmp_g2 = big_colstats(G2)


library(dplyr)
#save bigsnp object:
#a list with 3 elements
obj.bigsnp = list(genotypes = G2, # genotypes, FBM object
                  map = tibble(snp = 1:ncol(G2)), # map, i.e. SNP info
                  fam = tibble(FID = 1:nrow(G2))) # fam, i.e. info on individuals
#saving the bigsnp object
snp_save(obj.bigsnp)

#you can read the data with:
G3 = snp_attach("test_fbm.rds") # OBS: READ THE .rds(!!!) object
#genotypes are then stored in the .bk file.
