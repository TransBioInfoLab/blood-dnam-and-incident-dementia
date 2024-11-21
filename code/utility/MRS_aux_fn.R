library(survival)
library(janitor)
library(dplyr)
library(tibble)
library(ggplot2)
library(glmnet)
library(doParallel)
library(foreach)
# ==================================================================================
# Get MRS (Methylation Risk Score)
# ==================================================================================
get_MRS <- function(cpg.weight, betaset, scale = F){
  
  # Match CpGs between weight vector and beta matrix
  common.cpg <- intersect(
    names(cpg.weight), rownames(betaset)
  )
  
  # Extract overlapping CpGs from the beta matrix and weight vector
  beta.common <- betaset[common.cpg,]
  cpg.weight.common <- cpg.weight[common.cpg]
  
  # Optionally scale the beta values if scale is TRUE
  if(scale) {
    beta.common <- t(scale(t(beta.common)))
  }
  
  # Calculate MRS as a weighted sum of the beta values
  MRS <- (cpg.weight.common %*% beta.common)[1,]
  
  # Return the MRS and the list of matched CpGs
  return(
    list(
      MRS = MRS,
      common.cpg = common.cpg
    )
  )
}

# ==================================================================================
# Elastic Net Model Fitting
# ==================================================================================
fit_glmnet <- function(cpgs, beta, y, k = 10, a = NULL, seed = 1080,
                       fixed_df = NULL, scale = F, center = F,
                       lambda = NULL,
                       type = "regression", ...){
  
  # Prepare the design matrix from the beta values of the selected CpGs
  x <- t(beta[cpgs,])
  
  # Add additional fixed predictors to the design matrix if specified
  if(!is.null(fixed_df)) {
    pf <- c(rep(1, ncol(x)), rep(0, ncol(fixed_df)))
    x <- cbind(x, fixed_df)
  } else {
    pf <- rep(1, ncol(x))
  }
  
  # Remove rows with missing values from both x and y
  x <- x[!is.na(y),]
  y <- y[!is.na(y)]
  
  # Optionally scale and/or center the design matrix
  if(scale) {
    x[,1:length(cpgs)] <- scale(x[,1:length(cpgs)])
  }
  if(center) {
    x[,1:length(cpgs)] <- scale(x[,1:length(cpgs)], scale = F)
  }
  x <- as.matrix(x) # Ensure x is a matrix for glmnet
  
  # Determine the range of alpha values for cross-validation
  if(is.null(a)){
    a <- seq(0, 1, 0.1) # Test alpha values from 0 (Ridge) to 1 (Lasso)
    
    # Perform cross-validation over alpha values
    doParallel::registerDoParallel(4) # Use parallel processing for efficiency
    search <- foreach(i = a, .combine = rbind) %do% {
      set.seed(seed)
      foldid <- sample(1:k, size = length(y), replace = TRUE) # Create folds
      cv <- cv.glmnet(x = x, y = y, lambda = lambda,
                      foldid = foldid, alpha = i, penalty.factor = pf, ...)
      data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
    }
    print(search)
    
    # Identify the optimal alpha and lambda based on minimum CV error
    if(type == "survival"){
      cv <- search[which.max(search$cvm), ]
      cv_min <- cv$cvm
    } else {
      cv <- search[which.min(search$cvm), ]
      cv_min <- cv$cvm
    }
    l <- cv$lambda.1se
    alpha <- cv$alpha
  } else {
    # If alpha is fixed, perform cross-validation only for lambda
    set.seed(seed)
    foldid <- sample(1:k, size = length(y), replace = TRUE)
    cv <- cv.glmnet(x = x, y = y, lambda = lambda, 
                    foldid = foldid, alpha = a, penalty.factor = pf, ...)
    l <- cv$lambda.1se
    alpha <- a
    cv_min <- cv$cvm[cv$lambda == l]
  }
  
  # Fit the Elastic Net model using the optimal lambda and alpha
  set.seed(seed)
  mod <- glmnet(x = x, 
                  y = y, 
                  lambda = l, 
                  alpha = alpha,
                  penalty.factor = pf,
                  ...)
  
  # Extract model coefficients
  if(type == "survival"){
    coeff <- as.numeric(coef(mod))
    names(coeff) <- colnames(x)
  } else {
    coeff <- as.numeric(coef(mod))
    names(coeff) <- c("intersect", colnames(x))
  }
  
  # Return the fitted model and associated parameters
  mod <- list(
    coeff = coeff,
    cv_min = cv_min,
    lambda_min = l,
    a = alpha,
    mod = mod
  )
  
  return(mod)
}
