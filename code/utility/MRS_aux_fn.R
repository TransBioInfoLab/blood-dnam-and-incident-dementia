library(survival)
library(janitor)
library(dplyr)
library(tibble)
library(ggplot2)
library(glmnet)
library(doParallel)
library(foreach)
# ==================================================================================
# Get MRS
# ==================================================================================
get_MRS <- function(cpg.weight, betaset, scale = F){
  
  ## Match CpGs 
  common.cpg <- intersect(
    names(cpg.weight), rownames(betaset)
  )
  
  beta.common <- betaset[common.cpg,]
  cpg.weight.common <- cpg.weight[common.cpg]
  
  # ## Scaling CpGs
  if(scale) {
    beta.common <- t(scale(t(beta.common)))
  }
 
  MRS <- (cpg.weight.common %*% beta.common)[1,]
  
  return(
    list(
      MRS = MRS,
      common.cpg = common.cpg
    )
  )
}
# ==================================================================================
# Elastic Net
# ==================================================================================
fit_glmnet <- function(cpgs, beta, y, k = 10, a = NULL, seed = 1080,
                       fixed_df = NULL, scale = F, center = F,
                       lambda = NULL,
                       type = "regression",...){
  
  x <- t(beta[cpgs,])
  
  if(!is.null(fixed_df)) {
    pf <- c(rep(1, ncol(x)), rep(0, ncol(fixed_df)))
    x <- cbind(x, fixed_df)
  } else {
    pf <- rep(1, ncol(x))
  }
  
  # remove NA as glmnet cannot handle missing value
  x <- x[!is.na(y),]
  y <- y[!is.na(y)]
  
  if(scale) {
    x[,1:length(cpgs)] <- scale(x[,1:length(cpgs)])
  }
  if(center) {
    x[,1:length(cpgs)] <- scale(x[,1:length(cpgs)], scale = F)
  }
  x <- as.matrix(x)
  #y <- scale(y)
  
  if(is.null(a)){
    a <- seq(0, 1, 0.1)
    # cross-validation
    doParallel::registerDoParallel(4)
    search <- foreach(i = a, .combine = rbind) %do% {
      set.seed(seed)
      foldid <- sample(1:k,size=length(y),replace=TRUE)
      cv <- cv.glmnet(x = x, y = y, lambda = lambda,
                      foldid = foldid, alpha = i, penalty.factor = pf, ...)
      data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
    }
    print(search)
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
    set.seed(seed)
    foldid <- sample(1:k,size=length(y),replace=TRUE)
    cv <- cv.glmnet(x = x, y = y, lambda = lambda, 
                    foldid = foldid, alpha = a, penalty.factor = pf,...)
    l <- cv$lambda.1se
    alpha <- a
    
    cv_min <- cv$cvm[cv$lambda == l]
    
  }
  
  # Fit model
  set.seed(seed)
  ridge <- glmnet(x = x, 
                  y = y, 
                  lambda = l, 
                  alpha = alpha,
                  penalty.factor = pf,
                  ...)
  
  if(type == "survival"){
    coeff <- as.numeric(coef(ridge))
    names(coeff) <- colnames(x)
  } else {
    coeff <- as.numeric(coef(ridge))
    names(coeff) <- c("intersect", colnames(x))
  }
  
  
  mod <- list(
    coeff = coeff,
    cv_min = cv_min,
    lambda_min = l,
    a = alpha,
    mod = ridge
  )
  
  return(mod)
  
}
