#######################################################################################################
# =================================================================================================== #
# Single CpG assocation test
# =================================================================================================== #
#######################################################################################################
# Load package
library(tidyverse)
library(survival)
library(janitor)
# ===================================================================================================
# Single CpG testing: cox regression
# ===================================================================================================
# Function for cox regression analysis
get_cox_coef <- function(cpg, pheno, time_var, event_var, adjust_var = NULL, ...){

  # make survival formula
  fo <-  paste0("Surv(", time_var, ",", event_var, ") ~ cpg")
  if(!is.null(adjust_var)){
    fo <- paste0(fo, " + ", paste0(adjust_var, collapse = "+"))
  }
  formula <- as.formula(fo)
  
  data <- data.frame(cpg = cpg, pheno)
  
  # check warning
  w <- check_warning(
    formula, 
    data,
    ...
  )
  
  # fit cox regression
  cox_mod <- coxph(
    formula,
    data = data,
    ...
  )
  
  # get summary and coefficient
  cox_coef <- summary(cox_mod)$coefficients %>% data.frame()
  
  coef_df <- cox_coef[grepl("cpg", rownames(cox_coef)),] %>% 
    clean_names()
  
  coef_df$warning <- w
  
  return(coef_df)
}
# Wrapper function
cox_coef <- function(beta, pheno, time_var, event_var, adjust_var = NULL, 
                     parallel = T, 
                     scale = T,
                     fdr_method = "fdr", 
                     save = F, 
                     dir.save = NULL,
                     cores = 10,
                     prefix = "Framingham", ...){
  
  if(parallel) doParallel::registerDoParallel(cores)
  
  results <- plyr::adply(
    beta,
    .margins = 1,
    .fun = function(cg){
      if(scale) cg <- scale(cg)
      suppressWarnings({
        get_cox_coef(
          cpg = cg, 
          pheno = pheno,
          time_var = time_var,
          event_var = event_var,
          adjust_var = adjust_var,
          ...
        )
      })
      
    }, .id = "cpg", .parallel = parallel
  )

  colnames(results)[1] <- "cpg"
  # Add fdr
  results <- add_fdr(results, method = fdr_method)
  
  if(save){
    write_csv(
      results,
      file.path(dir.save, paste0(prefix, "_single_cpg_cox_results.csv"))
    )
  }
  return(results)
}

# Check results that exist warning
check_warning <- function(formula, data,...){
  
  mod <- tryCatch({
    coxph(
      formula,
      data = data,
      ...
    )
  }, warning = function(w) return(NULL))
  
  if(length(class(mod) > 2)){
    w <- F
  } else if (!is.null(mod) && class(mod) == "coxph"){
    w <- F
  } else {w <- T}
  
  return(w)
}

# Add fdr
add_fdr <- function(results, method = "fdr"){
  
  pr <- results[[grep("^p",colnames(results), value = T)]]
  results$fdr <- p.adjust(pr, method = method)
  
  return(results)
  
}
# ===================================================================================================
# Single CpG testing: covariates analysis
# ===================================================================================================
get_multi_coef <- function(cpg, pheno, test_var, adjust_var = NULL, ...){
  
  fo <- paste0("M ~ ", paste0(test_var, collapse = "+"))
  # make formula
  if(!is.null(adjust_var)){
    fo <- paste0(fo, "+", paste0(adjust_var, collapse = "+"))
  }
  formula <- as.formula(fo)
  
  data <- data.frame(M = minfi::logit2(cpg), pheno)
  
  # fit linear regression
  mod <- lm(
    formula,
    data = data,
    ...
  )
  ll <- levels(data[[test_var]])
  if(!is.null(ll) && length(ll) > 2) {
    fo <- paste0("M ~ ", paste0(adjust_var, collapse = "+"))
    formula <- as.formula(fo)
    mod0 <- lm(
      formula,
      data = data[!is.na(data[[test_var]]),],
      ...
    )
    
    mod_coef <- data.frame(
      p = anova(mod0, mod, test="Chisq")[2,5]
    )
    colnames(mod_coef) <- paste0(test_var, "_",colnames(mod_coef))
    
  } else {
    # get summary and coefficient
    mod_coef <- data.frame(summary(mod)$coefficients)
    mod_coef <- mod_coef[grep(paste0(test_var, collapse = "|"),rownames(mod_coef)),] %>% 
      janitor::clean_names() %>% 
      dplyr::select(estimate, pr_t) %>%
      dplyr::mutate(p = pr_t, .keep = "unused") %>%
      rownames_to_column("factors") %>%
      pivot_wider(names_from = c("factors"), 
                  values_from = !contains("factors"),
                  names_glue = "{factors}_{.value}",
                  names_vary = "slowest") 
  }
 
  return(mod_coef)
}
# ===================================================================================================
# Single CpG testing: partial correlation
# ===================================================================================================
partial_corr <- function(formula, data, dep_var = NULL, method = "spearman"){
  
  formula <- as.formula(formula)
  
  # split formula into variables
  vars <- all.vars(formula)
  
  # get independent variables
  res_var <- vars[1]
  if(is.null(dep_var)){
    # set the second variable as independent variable if dep_var is NULL
    dep_var <- vars[2] 
  }
  
  # select variables for test
  data_sel <- na.omit(data %>% dplyr::select(all_of(vars))) # remove NA
  
  # fit partial spearman correlation
  par_corr <- ppcor::pcor(data_sel, method = method)
  
  # Create summary 
  par_summ <- data.frame(
    var = vars[-1],
    rho = par_corr$estimate[!vars %in% res_var,res_var],
    statistic = par_corr$statistic[!vars %in% res_var,res_var],
    p.value = par_corr$p.value[!vars %in% res_var,res_var],
    n = par_corr$n
  )
  
  par_ind <- par_summ %>% filter(var == dep_var)
  par_ind$var <- NULL
  
  return(
    par_ind
  )
  
}

# ===================================================================================================
# Adjust M values function
# ===================================================================================================
methyl_adj <- function(mat, 
                       pheno, 
                       adjust_var, 
                       convert_to_M = T, 
                       return_to_beta = T,
                       parallel = T){
  
  if(convert_to_M){
    # transform beta to M
    mat <- minfi::logit2(mat)
  }
  mat <- as.matrix(mat)
  if(parallel) doParallel::registerDoParallel(10)
  
  resid_mat <- plyr::aaply(
    mat,
    .margins = 1,
    .fun = function(m){

      dat <- data.frame(M = m, pheno)
      # Create formula
      fo <- paste0("M ~ ", paste0(adjust_var, collapse = "+"))
      
      # Fit LM model
      lm_mod <- lm(as.formula(fo), data = dat)
      
      # Get residuals
      r <- resid(lm_mod)
      
      return(r)
    },.parallel = parallel
  )
  
  if(return_to_beta){
    # Transform back to beta
    resid_mat <- minfi::ilogit2(resid_mat)
  }
 
  return(resid_mat)
  
}
