library(mvtnorm)
library(MatchIt)
library(mice)



CCA_treatment_estimation <- function(res_row, res_col_name) {
  # Remove rows with NAN
  datCCA <- na.omit(dat)
  # Do matching
  m.out <- matchit(trt ~ X1+X2, data = datCCA, method = 'nearest',
                   ratio = 1) 
  dataMatched <- match.data(m.out)
  
  # Treatment effect
  resultsData[res_row, res_col_name] <<- mean(dataMatched$y[dataMatched$trt == 1]) - mean(dataMatched$y[dataMatched$trt == 0])
}

MI_treatment_estimation <- function(num_imputations, res_row, W_res_col_name, A_res_col_name) {
  m = num_imputations
  withinResults <- c()
  acrossResults <- c()
  propensityScores <- data.frame(matrix(ncol=0, nrow = n))
  if (imputation_index == 0){
    mids <- mice(dat[,c(3,4)], m = m, method = 'norm', print = FALSE)   
  }else if (imputation_index == 1){
    mids <- mice(dat[,c(1,3,4)], m = m, method = 'norm', print = FALSE)  
  }
  
  for(k in 1:m) {
    # WITHIN
    # Access imputed dataset
    dataImputed <- complete(mids, k)
    # Append y and trt
    dataImputed$trt = dat$trt
    dataImputed$y = dat$y
    # Do matching
    m.out <- matchit(trt ~ X1+X2, data = dataImputed, method = 'nearest',
                     ratio = 1)
    dataMatched <- match.data(m.out)
    # Append propensity scores of each obs in a dataset
    propensityScores <- cbind(propensityScores, data.frame(m.out$distance))
    
    # Treatment effect
    withinResults[k] <- mean(dataMatched$y[dataMatched$trt == 1]) - mean(dataMatched$y[dataMatched$trt == 0])
  }
  # ACROSS
  avg_prop <- rowMeans(propensityScores)
  m.out <- matchit(trt ~ X1+X2, data = dataImputed, method = 'nearest',ratio = 1, distance = avg_prop)
  dataMatched <- match.data(m.out)
  # Treatment effect
  acrossResults[1] <- mean(dataMatched$y[dataMatched$trt == 1]) - mean(dataMatched$y[dataMatched$trt == 0])
  
  # Store results
  resultsData[res_row, W_res_col_name] <<- mean(withinResults)  
  resultsData[res_row, A_res_col_name] <<- mean(acrossResults)  
}

#0 - MCAR
#1 - MAR1 - missingness only on control x2
#2 - MAR2 - missing on both trt and control x2
missingness_index <- 0

# 0 - assignemnt depends only on x1
# 1 - assignemnt depends only on x2
# 2 - assignment depends on both x1 and x2
treatment_index <- 0

#0 - do not include Y in the imputation
#1 - include Y in the imputation 
imputation_index <- 0

#0 - trt effects depends on X1
#1 - trt effects depends on X2
#2 - trt effects depends on X1 and X2
treatment_effect_index = 0

# Create a dataframe to store estimated treatment effects
resultsData <- data.frame(matrix(ncol = 0, nrow = 0))

res_std <- res <- data.frame(matrix(ncol = 12, nrow = 0))


n <- 1100
v <- 0
for (z in 0:2){
  missingness_index <- z
  
  for (l in 0:2){
    treatment_index <- l
      
    for (w in 0:2){
      treatment_effect_index = w
      
      for (r in 0:1){
        imputation_index <- r
        
          # start <- Sys.time()
          for (j in 1:250){
              print(v)
              # Generate covariates
              Sigma <- matrix(2.5, ncol = 2, nrow = 2)
              diag(Sigma) <- 5
              mu_X <- rep(10,2)
              
              X_raw <- X <- rmvnorm(n, mu_X, Sigma)
              colnames(X_raw) <- c('X1', 'X2')
              colnames(X) <- c('X1', 'X2')
              
              # Generate treatments
              #Trt assignment depends on x1
              if (treatment_index == 0){
                beta_trt <- 0.5
                beta0 <- -7.8
                xb <- beta0 +  beta_trt * X[,1] 
                p <- exp(xb) / (1 + exp(xb))
                trt <- rbinom(n, 1, p)  
              }
              
              #Trt assignment depends on x2
              if (treatment_index == 1){
                beta_trt <- 0.5
                beta0 <- -7.8
                xb <- beta0 +  beta_trt * X[,2] 
                p <- exp(xb) / (1 + exp(xb))
                trt <- rbinom(n, 1, p)  
              }
              
              #Trt assignment depends on x1 and x2
              if (treatment_index == 2){
                beta_trt <- 0.255
                beta0 <- -7.8
                xb <- beta0 + X %*% rep(beta_trt,2)
                p <- exp(xb) / (1 + exp(xb))
                trt <- rbinom(n, 1, p)  
              }
              
              # Generate responses 
              beta <- rep(1,2)
              beta0 <- 0
              
              #treatment effect only depends on x1
              if (treatment_index == 0){
                alpha0 <- 2.5
                alpha1 <- 0.25
                gamma <- alpha0 + alpha1*X[,1]
              }
              
              #treatment effect only depends on x2
              if (treatment_index == 1){
                alpha0 <- 2.5
                alpha1 <- 0.25
                gamma <- alpha0 + alpha1*X[,2]
              }
              
              #treatment effect depends on x1 and x2
              if (treatment_index == 2){
                alpha0 <- 2.5
                alpha1 <- 0.125
                alpha2 <- 0.125
                gamma <- alpha0 + alpha1*X[,1] + alpha2*X[,2]
              }
              
              mu1 <- beta0 + gamma + X %*% beta
              mu0 <- beta0 + X %*% beta
              sigma <- 1
              
              # Treatment
              y1 <- rnorm(n, mu1, sigma)
              # Control
              y0 <- rnorm(n, mu0, sigma)
              
              # Remove the unobserved treatments
              y <- y0
              y[trt == 1] <- y1[trt == 1]
              
              #Check baseline bias
              #True effect
              delta <- mean(y1) - mean(y0)
              resultsData[j, 'Complete.Data'] = delta     
              
              #Biased effect.  
              mean(y1[trt == 1]) - mean(y0[trt == 0])
              
              #Add missingness in X1
              #MCAR
              if (missingness_index == 0) {
                ind <- sample(1:nrow(X), nrow(X) / 2, replace = FALSE)
                X[ind, 2] <- NA
              }
              
              #MAR1
              if (missingness_index == 1) {
                logi <- -10.1 + 0.9 * X[, 1]
                p <- exp(logi) / (1 + exp(logi))
                R <- rbinom(nrow(X), 1, p)
                X[R == 1 & trt ==0, 2] <- NA
              }
              
              #MAR2
              if (missingness_index == 2) {
                logi <- -10.1 + 0.9 * X[, 1]
                p <- exp(logi) / (1 + exp(logi))
                R <- rbinom(nrow(X), 1, p)
                X[R == 1 & trt==0, 2] <- NA
                
                X[trt == 1 & rbinom(n,1,0.3), 2] <- NA
              }
      
              # Create a dataframe with missingness
              dat <- data.frame(y, trt, X)
              
              #CCA
              CCA_treatment_estimation(j, 'CCA')
              #MI 5, 10, 15, 20, 50
              MI_treatment_estimation(5, j,  'W_MI5', 'A_MI5')
              MI_treatment_estimation(20, j, 'W_MI20', 'A_MI20')
              MI_treatment_estimation(50, j, 'W_MI50', 'A_MI50')
          }
          # end <- Sys.time()
        results <- t(data.frame(colMeans(resultsData)))
        temp <- data.frame(Missingness_Mech = z, Trt_Assignment = l, Trt_Effect_Size = w, Include_Y_in_Imputation = r, results)
        res <- rbind(res, temp)
        
        stdev <- t(sqrt(apply(resultsData, 2, var)))
        temp <- data.frame(Missingness_Mech = z, Trt_Assignment = l, Trt_Effect_Size = w, Include_Y_in_Imputation = r, stdev)
        res_std <- rbind(res_std, temp)
        v <- v +1

      }
    }
  }
}