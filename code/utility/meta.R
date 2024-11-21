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
# meta analysis with metagen
# ===================================================================================================
# ===================================================================================================
# Meta-analysis Wrapper Function
# ===================================================================================================
# This function performs a meta-analysis on a list of results from multiple studies.
meta_wrapper <- function(results_list, effect, se, full_table = T, sm = "HR", test_var = "cpg",
                         select_var = c("seqnames", "start", "end", "width", "Islands.UCSC.Relation_to_Island",
                                        "UCSC_RefGene_Name", "UCSC_RefGene_Group", "GREAT_annotation")) {
  
  # Identify CpGs present in at least two studies
  cpg <- unlist(purrr::map(results_list, ~ .[[test_var]]))
  cpg_tb <- table(cpg)
  cpg_select <- names(cpg_tb[cpg_tb > 1])
  
  # Register parallel backend for improved performance
  doParallel::registerDoParallel(30)
  
  # Match results to the selected common CpGs and perform meta-analysis
  results <- foreach::foreach(
    i = cpg_select,
    .errorhandling = "remove"
  ) %dopar% {
    dat <- purrr::map(results_list, ~ .[.[[test_var]] == i, ])
    dat <- purrr::reduce(dat, rbind)
    dat <- purrr::compact(dat)
    
    # Perform meta-analysis on the aggregated data
    meta_fn(dat = dat, 
            effect = effect,
            se = se,
            full_table = full_table,
            sm = sm,
            test_var = test_var,
            select_var = select_var)
  }
  
  # Stop the parallel backend
  doParallel::stopImplicitCluster()
  
  # Combine results from all parallel tasks into a single data.table
  results <- data.table::rbindlist(results, fill = T)
  
  # Add FDR-corrected p-values
  results <- mutate(results, 
                    fdr = p.adjust(pVal.fixed, method = 'fdr'), .before = k)
  
  return(results)
}

# ===================================================================================================
# Meta-analysis Function
# ===================================================================================================
# This function performs meta-analysis on a dataset for a single CpG.
meta_fn <- function(dat, effect, se, full_table = T, sm = "HR", return_metagen = F, test_var = "cpg", 
                    select_var = c("seqnames", "start", "end", "width", "Islands.UCSC.Relation_to_Island",
                                   "UCSC_RefGene_Name", "UCSC_RefGene_Group", "GREAT_annotation"),
                    ...) {
  
  # Perform the meta-analysis using the metagen function from the meta package
  f <- metagen(
    TE = dat[[effect]],
    seTE = dat[[se]],
    sm = sm,
    ...
  )
  
  # Compile the fixed- and random-effect meta-analysis results into a data frame
  mod_coef <- data.frame(
    cpg = unique(dat[[test_var]]),
    estimate = f$TE.fixed,
    se = f$seTE.fixed,
    pVal.fixed = f$pval.fixed,
    pVal.random = f$pval.random,
    pVal.Q = f$pval.Q,
    direction = paste0(ifelse(dat[[effect]] > 0, "+", "-"), collapse = ""),
    k = nrow(dat)
  )
  colnames(mod_coef)[1] <- test_var
  
  # If full_table is TRUE, include additional columns from the input data
  if(full_table) {
    dat$warning <- NULL
    dat$study <- NULL
    col_names <- c(test_var, select_var)
    dat_wider <- pivot_wider(
      dat, 
      id_cols = all_of(col_names),
      names_from = study_id,
      names_glue = "{study_id}_{.value}",
      values_from = colnames(dat)[!colnames(dat) %in% c(col_names, "study_id")]
    )
    mod_coef <- left_join(mod_coef, dat_wider, by = test_var)
  } 
  
  # Optionally return the full metagen object instead of the summary table
  if(return_metagen) {
    return(f)
  } else {
    return(mod_coef)
  }
}
