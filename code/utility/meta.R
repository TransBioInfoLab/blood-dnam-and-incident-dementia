#######################################################################################################
# ===================================================================================================
# Function for meta analysis 
# ===================================================================================================
#######################################################################################################
library(metafor)
library(meta)
library(survival)
library(doParallel)
# ===================================================================================================
# meta analysis with metafor
# ===================================================================================================
metafor_wrapper <- function(results_list, rho_df, effect, se, full_table = T, robust = T){
  
  # Get cpgs that have at least two study
  cpg <- unlist(purrr::map(results_list, ~.[["cpg"]]))
  cpg_tb <- table(cpg)
  cpg_select <- names(cpg_tb[cpg_tb > 1])
  
  doParallel::registerDoParallel(30)
  # Match to the common cpg
  results <- foreach::foreach(
    i = cpg_select,
    .errorhandling = "remove"
  ) %dopar% {
    dat <- purrr::map(results_list, ~.[.$cpg == i,])
    dat <- purrr::reduce(dat, rbind) 
    dat <- purrr::compact(dat)

    if(is.null(nrow(rho_df))) {
      rho <- rho_df
    } else {
      if(unique(dat$cpg) %in% rho_df$cpg) {
        rho <- rho_df[rho_df$cpg %in% unique(dat$cpg),"corr"]
      }
    }
   
    # Meta analysis
    metafor_fn(dat = dat, 
               rho = rho,
               effect = effect,
               se = se,
               full_table = full_table,
               robust = robust)
  }
  
  doParallel::stopImplicitCluster()
  
  results <- data.table::rbindlist(results, fill = T)
  results <- mutate(results,
                    fdr = p.adjust(pval, method = 'fdr'), .before = k)
  results
  
}

metafor_fn <- function(dat, rho = NA, effect, se, full_table = T, robust = T){
  
  study_num <- table(dat$study)
  
  if(study_num["FHS"] == 2 && !is.na(study_num["FHS"])) {
    V <- vcalc(vi = dat[[se]]^2, cluster = study, obs = study_id, rho = rho, data = dat)
  } else {
    V <- vcalc(vi = dat[[se]]^2, cluster = study, data = dat)
    rho <- NA
  }
  
  mod <- rma.mv(yi = dat[[effect]], V = V, data = dat)
  if(robust) {
    mod <- robust(mod, cluster = study, clubSandwich = T)
  }
 
  mod_coef <- data.frame(cpg = unique(dat$cpg),
                         coefficients(summary(mod)),
                         rho = rho,
                         direction = paste0(ifelse(dat[[effect]] > 0, "+", "-"), 
                                             collapse = ""),
                         k = nrow(dat))
  
  if(full_table) {
    # column to retain
    dat$warning <- NULL
    dat$study <- NULL
    col_names <- c("cpg", "seqnames", "start", "end", "width", "Islands.UCSC.Relation_to_Island",
                   "UCSC_RefGene_Name", "UCSC_RefGene_Group", "GREAT_annotation")
    dat_wider <- pivot_wider(
      dat, id_cols = all_of(col_names),
      names_from = study_id,
      names_glue = "{study_id}_{.value}",
      values_from = colnames(dat)[! colnames(dat) %in% c(col_names, "study_id")]
    )
    
    mod_coef <- left_join(mod_coef, dat_wider, by = "cpg")
  } 

  mod_coef
  
}

calc_corr <- function(beta_list, common_cpg, method = "spearman"){
  
  doParallel::registerDoParallel(10)

  results <- plyr::laply(
    common_cpg,
    .fun = function(i){
    dat <- purrr::map(beta_list, ~.[i,])

    # correlation
    cor(
      dat[[1]], dat[[2]],
      method = method
    )
    
  })
  
  doParallel::stopImplicitCluster()
  names(results) <- common_cpg
  results

  
}
# ===================================================================================================
# meta analysis with metagen
# ===================================================================================================
meta_wrapper <- function(results_list, effect, se, full_table = T, sm = "HR", test_var = "cpg",
                         select_var = c("seqnames", "start", "end", "width", "Islands.UCSC.Relation_to_Island",
                                        "UCSC_RefGene_Name", "UCSC_RefGene_Group", "GREAT_annotation")){
  
  # Get cpgs that have at least two study
  cpg <- unlist(purrr::map(results_list, ~.[[test_var]]))
  cpg_tb <- table(cpg)
  cpg_select <- names(cpg_tb[cpg_tb > 1])
  
  doParallel::registerDoParallel(30)
  # Match to the common cpg
  results <- foreach::foreach(
    i = cpg_select,
    .errorhandling = "remove"
  ) %dopar% {
    dat <- purrr::map(results_list, ~.[.[[test_var]] == i,])
    dat <- purrr::reduce(dat, rbind) 
    dat <- purrr::compact(dat)
    
    # Meta analysis
    meta_fn(dat = dat, 
            effect = effect,
            se = se,
            full_table = full_table,
            sm = sm,
            test_var = test_var,
            select_var = select_var)
  }
  
  doParallel::stopImplicitCluster()
  
  results <- data.table::rbindlist(results, fill = T)
  results <- mutate(results, 
                    fdr = p.adjust(pVal.fixed, method = 'fdr'), .before = k)
  results
  
}

meta_fn <- function(dat, effect, se, full_table = T, sm = "HR", return_metagen = F, test_var = "cpg", 
                    select_var = c("seqnames", "start", "end", "width", "Islands.UCSC.Relation_to_Island",
                                   "UCSC_RefGene_Name", "UCSC_RefGene_Group", "GREAT_annotation"),
                    ...){
  
  f <- metagen(
    TE = dat[[effect]],
    seTE = dat[[se]],
    sm = sm,
    ...
  )
  
  mod_coef <- data.frame(
    cpg = unique(dat[[test_var]]),
    estimate = f$TE.fixed,
    se = f$seTE.fixed,
    pVal.fixed = f$pval.fixed,
    pVal.random = f$pval.random,
    pVal.Q = f$pval.Q,
    direction = paste0(ifelse(dat[[effect]] > 0, "+", "-"), 
                       collapse = ""),
    k = nrow(dat)
  )
  colnames(mod_coef)[1] <- test_var
  if(full_table) {
    # column to retain
    dat$warning <- NULL
    dat$study <- NULL
    col_names <- c(test_var, select_var)
    dat_wider <- pivot_wider(
      dat, id_cols = all_of(col_names),
      names_from = study_id,
      names_glue = "{study_id}_{.value}",
      values_from = colnames(dat)[! colnames(dat) %in% c(col_names, "study_id")]
    )
    
    mod_coef <- left_join(mod_coef, dat_wider, by = test_var)
  } 
  
  if(return_metagen) {
    return(f)
  } else {
    return(mod_coef)
  }
  
}