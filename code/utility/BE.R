####### Run BEclear to remove batch effects #######
library(BEclear)
rm_batch_BE <- function(se, BE_summ = NULL, BE_score = NULL, batch_var, sample_var, cores = 10, prefix = "FHS_EPIC", dir.data, dir.results, ...){
  
  # Retrieve betaMat
  betaMat <- assay(se)
  
  # Retrieve pheno data
  pheno_df <- colData(se) %>% 
    as.data.frame() %>%
    mutate(batch_id = as.factor(get(batch_var)),
           sample_id = get(sample_var))

  if(is.null(BE_summ)){
    message("Start pvalue calculation")
    # Calculate batch pvalue
    doParallel::registerDoParallel(cores)
    BE <- calcBatchEffects(
      data = betaMat,
      samples = pheno_df,
      BPPARAM = MulticoreParam(cores)
    )
    
    # Save results
    save(
      BE,
      file = file.path(paste0(dir.data, "/", prefix, "_BE.rda"))
    )
    
    # Summarize median comparison
    BE_summ <- calcSummary(BE$med, BE$pval)
    
    # Calculate score
    BE_score <- calcScore(
      data = betaMat,
      samples = pheno_df,
      BE_summ
    )
    
    # Save results
    writexl::write_xlsx(
      list(summary = BE_summ, score = BE_score),
      file.path(paste0(dir.results, "/", prefix, "_summ_scores.xlsx"))
    )
  }
  
  # Remove probes in batch with BEscores > 0.02 
  rm_batch <- BE_score[BE_score$BEscore > 0.02, ]
  rm_batch <- rm_batch$batch_id
  
  if(length(rm_batch) > 0) {

    message("Number of Batch > 0.02", length(rm_batch))
    summ <- BE_summ[BE_summ$batch_id %in% rm_batch,]
    clearedMatrix <- clearBEgenes(betaMat, pheno_df, summ)
    
    predicted <- imputeMissingData(
      data = clearedMatrix,
      BPPARAM = MulticoreParam(cores),
      ...
    )
    
    save(
      predicted,
      file = file.path(paste0(dir.data, "/", prefix, "_BE_predicted.rda"))
    )
  }
  
}