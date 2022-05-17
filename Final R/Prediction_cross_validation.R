#' Description - Function to populate the FBM
#' 
#' @param person
#' @param k
#' @param threshold
#' @param disease
#' @param method
#' @param liabilities
#' @example 
#' 

Prediction_cross_validation <- function(person, k = 10, threshold = 0, disease, method = "GWAS", liabilities = person$FAM$Status) {
  n <- nrow(person$genotypes)
  #res <- tibble(threshold = numeric(), PRS = vector(), predictive_power = numeric())
  block_size <- n%/%k
  bestest_score <- 0
  bestest_model <- NULL
  best_scores <- numeric(length(threshold))
  avg_scores <- numeric(length(threshold))
  for (i in 1:length(threshold)){
    #print(i)
    scores <- numeric(k)
    best_score <- 0
    for (j in 0:(k-1)){
      #print(j)
      
      #Create blocks
      block_start <- j*block_size + 1
      if (j == k-1){
        block_end <- n} else {
      block_end <- (j+1)*block_size}
      
      #Train on k-1 folds
      if(method == "GWAS"){
        regr <- GWAS(person, person$FAM$Status, include = c(1:n)[-(block_start:block_end)])
      }
      else if(method == "GWAX"){
        regr <- GWAX(person, include = c(1:n)[-(block_start:block_end)])
      }
      else if(method == "LTFH"){
        regr <- GWAS(person, liabilities, include = c(1:n)[-(block_start:block_end)])
      }
      
      #Calculate PRS on 1 fold
      PRS <- snp_PRS(G = person$genotypes, betas.keep = regr$estim, ind.test = block_start:block_end, lpS.keep = -log10(regr$p.value), thr.list = threshold[i])
      
      #Normalize PRS with mean = 0 and sd = 1
      normalized_PRS <- (PRS - mean(PRS))/sd(PRS)
      
      #Find score of model
      score <- cor(normalized_PRS, person$FAM$Full_Liability[block_start:block_end])
      scores[j+1] <- score
      if(score > best_score){
        best_score <- score
        best_model <- regr
      }
      #print(scores)
      
    }
    #Update best models, scores
    best_scores[i] <- best_score[1,1]
    avg_scores[i] <- mean(scores)
    if(best_score > bestest_score){
      bestest_score <- best_score[1,1]
      bestest_model <- best_model
    }
  }
  results <- tibble(Pvalue = threshold, Average_Score = avg_scores, Best_Score = best_scores)
  return(list(Results = results, Best_Score = bestest_score, Best_Model = data.frame(bestest_model)))
}

#reg <- lm(person$FAM$Status~normalized_PRS)

#nn <- c(1,2,3,4,5)
#oo <- c(5,4,3,2,6,6,42,5,6,7)
#op <- list(c(5,4,3,2,6,6,42,5,6,7))
#print(op)

#mm <- 0.95

#a <- tibble(pvals = nn, PRS = list(vector(length = 10)), predpower = NA,)
#a
#a[1,3] <- mm
#a[1,2] <- op
#a


#b <- tibble(pvals = numeric(), PRS = vector(), predpower = numeric())
#b

#c <- as_tibble_row(list(pvals = nn[1], PRS = op, predpower = mm))

#rbind(b, c) 

#print(list(rnorm(10)))
#add_row(b, list(pvals = nn[1], PRS = list(rnorm(10)), predpower = mm))


#test <- data.frame("a" = z, "b" = aa)
#test[1,] <- c(5, c(1,2,3))
# <- as.list(5, c(1,2,3))
#est
#z <- c(0,1,1,0,1)
#aa <- c(1,0,1,1,1)
#sum(z == aa)
#Prediction_cross_validation(persona, k=3)
