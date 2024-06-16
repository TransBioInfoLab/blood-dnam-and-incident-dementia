library(survival)
library(janitor)
library(dplyr)
library(tibble)
library(pROC)
library(ggplot2)
# ==================================================================================
# Data preparation
# ==================================================================================
data_prep <- function(Data = "Framingham", dir.data, remove_65 = F){
  
  dir.data.selected <- file.path(dir.data, Data) 
  dir.data.pca <- file.path(dir.data.selected, "DNAm/pca_filtering/") 
  
  if(Data == "AIBL"){
    
    aibl.se <- readRDS(
      file.path(dir.data.pca, "AIBL_QNBMIQ_PCfiltered_age_at_least_65.RDS")
    )
    
    aibl.se <- aibl.se[,!aibl.se$`disease status:ch1` %in% "Mild Cognitive Impairment" ]
    table(aibl.se$`disease status:ch1`)
    # Alzheimer's disease     healthy control 
    #  135                    356 
    pheno <- colData(aibl.se)
    pheno$DIAGNOSIS <- factor(pheno$`disease status:ch1`, levels = c("healthy control", "Mild Cognitive Impairment", "Alzheimer's disease"))
    pheno$SEX <- ifelse(pheno$Sex == "M","Male","Female") %>% factor(., levels = c("Female", "Male"))
    pheno$AGE <- pheno$age.pred.Elastic_Net
    pheno$DIAGNOSIS <- pheno$DIAGNOSIS %>% droplevels()
    
    beta = assay(aibl.se)
  }
  
  if(Data == "ADNI"){
    
    adni.se <- readRDS(
      file.path(dir.data.pca, "ADNI_QNBMIQ_PCfiltered_min_age_at_visit_65_with_smoke_probes.RDS")
    )
    adni.se <- adni.se[,!adni.se$DX %in% "MCI" ]
    adni.se <- adni.se[,!is.na(adni.se$DX)]
    
    pheno <- colData(adni.se)
    pheno$DIAGNOSIS <- factor(
      ifelse(pheno$DX == "CN", "healthy control", "Alzheimer's disease"),
      levels = c("healthy control", "mild cognitive impairment", "Alzheimer's disease")
    ) %>% droplevels()
    pheno$SEX <- pheno$PTGENDER %>% 
      factor(., levels = c("Female", "Male"))
    pheno$AGE <- pheno$age_at_visit
    
    # Keep only last visit 
    pheno <- pheno %>% 
      as.data.frame %>%  
      dplyr::group_by(RID) %>% 
      filter(age_at_visit == max(age_at_visit)) 
    
    rownames(pheno) <- pheno$barcodes
    
    adni.se <- adni.se[,adni.se$barcodes %in% pheno$barcodes]
    adni.se <- adni.se[,pheno$barcodes]
    
    beta <- assay(adni.se)
    
  }
  
  if(Data == "AddNeuroMed"){
    
    addNeuroMed.se <- readRDS(file.path(dir.data.pca, "addNeuroMed_QNBMIQ_PCfiltered.RDS"))
    
    addNeuroMed.se <- addNeuroMed.se[,addNeuroMed.se$disease.state.ch1 !=  "mild cognitive impairment"]
    addNeuroMed.se <- addNeuroMed.se[rowSums(is.na(assay(addNeuroMed.se))) != ncol(addNeuroMed.se),]
    
    #  0 = control, 1 = AD
    pheno <- data.frame(colData (addNeuroMed.se)) 
    pheno$DIAGNOSIS <- factor(pheno$disease.state.ch1, levels = c("control", "mild cognitive impairment", "Alzheimer's disease"))
    pheno$SEX <- pheno$Sex.ch1 %>% 
      factor(., levels = c("Female", "Male"))
    pheno$DIAGNOSIS <- pheno$DIAGNOSIS %>% droplevels()
    pheno$AGE <- pheno$age.ch1
    plyr::count(pheno$DIAGNOSIS)
    
    addNeuroMed.se <- addNeuroMed.se[,pheno$geo_accession]
    
    beta <- assay(addNeuroMed.se)
  }
  
  if(Data == "Framingham"){
    
    framingham.se <- readRDS(file.path(dir.data.pca, "Framingham_Offspring_Exam8_se_exclude_DEM_JHU_UMN.RDS"))
    pheno <- data.frame(colData (framingham.se)) 
    pheno$DIAGNOSIS <- factor(ifelse(pheno$DEM_STATUS == 1, "DEM", "CN"), levels = c("CN", "DEM"))
    pheno$AGE <- pheno$AGE_AT_VISIT
    pheno$SEX <- ifelse(pheno$SEX == "female", "Female", "Male") %>% 
      factor(., levels = c("Female", "Male"))
    
    pheno <- plyr::adply(
      pheno,
      .margins = 1, 
      .fun = function(cl){
        cl$DEM_FOLLOW_UP_AT_Y5 <- ifelse(cl$DEM_FOLLOW_UP_Y_SINCE_EXAM8 <= 5 & cl$DEM_STATUS == 1, "DEM", "CN")
        cl$DEM_FOLLOW_UP_AT_Y10 <- ifelse(cl$DEM_FOLLOW_UP_Y_SINCE_EXAM8 <= 10 & cl$DEM_STATUS == 1, "DEM", "CN")
        cl$DEM_FOLLOW_UP_AT_Y14 <- ifelse(cl$DEM_FOLLOW_UP_Y_SINCE_EXAM8 <= 14 & cl$DEM_STATUS == 1, "DEM", "CN")
        cl$AD_FOLLOW_UP_AT_Y5 <- ifelse(cl$DEM_FOLLOW_UP_Y_SINCE_EXAM8 <= 5 & cl$AD_STATUS == 1, "DEM", "CN")
        cl$AD_FOLLOW_UP_AT_Y10 <- ifelse(cl$DEM_FOLLOW_UP_Y_SINCE_EXAM8 <= 10 & cl$AD_STATUS == 1, "DEM", "CN")
        cl$AD_FOLLOW_UP_AT_Y14 <- ifelse(cl$DEM_FOLLOW_UP_Y_SINCE_EXAM8 <= 14 & cl$AD_STATUS == 1, "DEM", "CN")
        return(cl)
      }
    )
    
    pheno <- pheno %>% 
      mutate(Age_range = ifelse(AGE < 65, "Age < 65", 
                                ifelse(AGE >= 65 & AGE < 75, "65 <= Age < 75", 
                                       ifelse(AGE >= 75 & AGE < 85, "75 <= Age < 85", "Age >= 85")))) %>%
      filter(DEM_FOLLOW_UP_Y_SINCE_EXAM8 >= 0)
    if(remove_65){
      pheno <- pheno %>% 
      filter(Age_range != "Age < 65") 
    }
    rownames(pheno) <- pheno$LABID
    beta <- assay(framingham.se)[,pheno$LABID]
    
  }
  
  return(list(
    beta = beta,
    pheno = pheno
  ))
  
}
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
# Test association
# ==================================================================================
test_assoc_mrs <- function(data_ls, 
                           model.formula = "DIAGNOSIS ~ MRS + AGE + SEX + B + NK + CD4T + Mono + Gran",
                           cpg.weight, 
                           DEM_var, 
                           source, 
                           scale = F,
                           test.datasets = "AddNeuroMed") {
  
  pheno <- data_ls[[test.datasets]]$pheno
  data <- data_ls[[test.datasets]]$beta
  
  MRS_ls <- get_MRS(cpg.weight, data, scale = scale)

  # print("Calculating MRS")
  # Calculate the MRS for each sample beta dot product weights
  MRS <- MRS_ls$MRS
  
  pheno$MRS <- MRS
  pheno$Gran <- pheno$Neutro + pheno$Eosino
  dem <- pheno[[DEM_var]]
  
  nsamples.testing <- nrow(pheno)
  common.cpgs <- MRS_ls$common.cpg
  
  model <- glm(
    model.formula,
    data = pheno, 
    family = binomial
  )

  model_summ <- coefficients(summary(model)) %>% 
    janitor::clean_names()
  pval <- model_summ[grep("mrs", rownames(model_summ)), "pr_z"]
  if(length(pval) == 0) pval <- NA
  
  probabilities <- predict(model, type = "response")

  pheno$probabilities <- probabilities
  
  res.roc <- roc(dem, probabilities)
  
  wilcox <- tryCatch({
    wilcox.test(as.formula(paste0("probabilities ~ ", DEM_var)), 
                data = pheno, exact = F, conf.int = TRUE)
  }, error = function(e){
    data.frame("p.value" = NA,"estimate" = NA)
  })
  
  ci <- tryCatch({
    ci.auc (res.roc, conf.level = 0.95, method = "bootstrap")
  }, error = function(e){
    NA
  })
  
  cutpoint <- coords (
    res.roc, 
    x = "best", 
    best.method = "youden",
    ret = c(
      "threshold", 
      "specificity", "sensitivity",
      "accuracy", "ppv", "npv"
    )
  )
  if(nrow(cutpoint) > 1) cutpoint <- cutpoint[1,]
  
  results <- cbind(
    data.frame(
      "Testing dataset" = test.datasets,
      "nsamples.testing" = nsamples.testing,
      "Model" = model.formula,
      "Number.Cpgs" = length(common.cpgs),
      "Cpgs" = source,
      "AUC" = res.roc$auc,
      "AUC CI" = paste0("95% CI: ",round(min(ci),digits = 4),"-" ,round(max(ci),digits = 4)," (DeLong)"),
      "Wilcoxon_pvalue_prob" = wilcox$p.value,
      "Wilcoxon_pvalue_estimate" = wilcox$estimate,
      "MRS_pval" = pval
    ), cutpoint
  )
  
  list(results = results,
       MRS = data.frame(probabilities = probabilities, MRS = MRS, DEM = dem))
  
}

