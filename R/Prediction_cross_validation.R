#' Perform Prediction Cross-validation to determine best model
#' 
#' This function is used to calculate predictive powers of different models.
#' 
#' @param rds.obj A .rds file with an FBM.code256 and accompanying FAM and MAP.
#' @param k Number of folds to be used in cross-validation. Default is 10.
#' @param threshold Vector of P-values to be used in thresholding. Default does not use thresholding.
#' @param disease List with properties of disease.
#' @param method Method to use for prediction. Possible methods are "GWAS", "GWAX", "LTFH". Default is "GWAS".
#' @param liabilities Vector of liabilities used for prediction with "LTFH" method. If not specified, uses "GWAS" method instead.
#' @return A list with 3 values: a tibble with average and best scores for each threshold, the best score, and the best model.
#' @export
#' 


Prediction_cross_validation <- function(rds.obj, k = 10, threshold = 0, disease, method = "GWAS", liabilities = person$FAM$Status) {
  n <- nrow(person$genotypes)
  block_size <- n%/%k
  bestest_score <- 0
  bestest_model <- NULL
  bestest_score1 <- 0
  bestest_score2 <- 0
  best_scores <- numeric(length(threshold))
  avg_scores <- numeric(length(threshold))
  #best_score1s <- numeric(length(threshold))
  #best_score2s <- numeric(length(threshold))
  #if(method == "GWAX"){
  #  proxy_status <- GWAX(person)$Proxy_Status
  #}
  
  for (i in 1:length(threshold)){
    scores <- numeric(k)
    best_score <- 0
    #best_score1 <- 0
    #best_score2 <- 0
    for (j in 0:(k-1)){
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
        regr <- GWAX(person, include = c(1:n)[-(block_start:block_end)])$GWAS_Data
      }
      else if(method == "LTFH"){
        regr <- GWAS(person, liabilities, include = c(1:n)[-(block_start:block_end)])
      }
      
      #Calculate PRS on 1 fold
      PRS <- bigsnpr::snp_PRS(G = person$genotypes, betas.keep = regr$estim, ind.test = block_start:block_end, lpS.keep = -log10(regr$p.value), thr.list = threshold[i])
      
      #Normalize PRS with mean = 0 and sd = 1
      normalized_PRS <- (PRS - mean(PRS))/sd(PRS)
      
      #Find score of model
      #Correlation between predicted fitted values and status = R2
      # if(method == "GWAS"){
      #   regr2 <- lm(person$FAM$Status[block_start:block_end] ~ normalized_PRS - 1)
      #   score1 <- cor(regr2$fitted.values, person$FAM$Status[block_start:block_end])^2
      #   print(cor(regr2$fitted.values, normalized_PRS))
      # }
      # else if(method == "GWAX"){
      #   regr2 <- lm(proxy_status[block_start:block_end] ~ normalized_PRS - 1)
      #   score1 <- cor(regr2$fitted.values, proxy_status[block_start:block_end])^2
      # }
      # else if(method == "LTFH"){
      #   regr2 <- lm(liabilities[block_start:block_end] ~ normalized_PRS - 1)
      #   score1 <- cor(regr2$fitted.values, proxy_status[block_start:block_end])^2
      # }
      #
      
      #Correlation between PRS and status
      score <- cor(normalized_PRS, person$FAM$Status[block_start:block_end])
      
      #score <- score1 + score2
      
      #Update best models, scores
      scores[j+1] <- score
      if(score > best_score){
        best_score <- score
        best_model <- regr
        #best_score1 <- score1
        #best_score2 <- score2
      }
    }
    #Update best models, scores
    best_scores[i] <- best_score#[1,1]
    avg_scores[i] <- mean(scores)
    #best_score1s[i] <- best_score1
    #best_score2s[i] <- best_score2
    if(best_score > bestest_score){
      bestest_score <- best_score#[1,1]
      bestest_model <- best_model
      #bestest_score1 <- best_score1
      #bestest_score2 <- best_score2
    }
  }
  results <- tibble::tibble(Pvalue = threshold, Average_Score = avg_scores, Best_Score = best_scores, R2 = best_scores^2)
  return(list(Results = results, Best_Score = bestest_score, Best_Model = data.frame(bestest_model)))
}