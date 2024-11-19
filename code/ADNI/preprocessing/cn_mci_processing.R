cohort <- "ADNI"
dir.base <- "."
dir.data <- file.path(dir.base,"../DATASETS/", cohort)
dir.data.pca <- file.path(dir.data, "DNAm/pca_filtering/all_sample") 
dir.data.pheno.raw <- file.path(dir.data, "Phenotype/raw/")
dir.data.clinical <- file.path(dir.data, "Phenotype/processed/all_sample")


adni.se <- readRDS(file.path(dir.data.pca, "SummarizedExperiment_PCfiltered.RDS"))

# Survival 
pheno <- colData(adni.se) %>% as.data.frame()
pheno_baseline <- pheno %>%
  group_by(RID) %>% 
  slice_min(age_at_visit)

# select CN and MCI samples at baseline
pheno_baseline_cn <- pheno_baseline %>%
  filter(DX %in% c("CN", "MCI"))
rid <- pheno_baseline_cn$RID
length(rid) # 543

# Read all clinical data and get the last diagnosis results
clinical <- readr::read_csv(file.path(dir.data.pheno.raw,"ADNIMERGE_downloaded_3-1-2021.csv"))
clinical_filter <- clinical %>%
  dplyr::select(RID, VISCODE, COLPROT,EXAMDATE, DX)
clinical_filter <- clinical_filter[!is.na(clinical_filter$DX),] # rm NA
clinical_filter <- clinical_filter %>%
  filter(RID %in% rid) %>%
  group_by(RID) %>%
  slice_max(EXAMDATE)

# Read latest diagnosis file and get the latest diagnosis results
diagnosis <- read_csv(
  file.path(dir.data.pheno.raw, "DXSUM_PDXCONV_ADNIALL_downloaded_12-11-2022.csv")
)
diagnosis_filter <- diagnosis %>%
  filter(RID %in% rid) 
diagnosis_filter$EXAMDATE <- as.Date(diagnosis_filter$EXAMDATE, "%Y/%m/%d") 
diagnosis_filter <- diagnosis_filter %>%
  group_by(RID) %>%
  slice_max(EXAMDATE) %>% 
  mutate(DX2 = ifelse(DIAGNOSIS == 1, "CN", ifelse(DIAGNOSIS == 2, "MCI", ifelse(DIAGNOSIS == 3, "Dementia", NA))),
         .keep = "unused")
diagnosis_filter <- diagnosis_filter[match(clinical_filter$RID, diagnosis_filter$RID),]

# update the last diagnosis results to latest diagnosis results
clinical_filter$DX2 <- clinical_filter$DX
clinical_filter$DX2[!is.na(diagnosis_filter$DX2)] <- diagnosis_filter$DX2[!is.na(diagnosis_filter$DX2)]
clinical_filter$EXAMDATE2 <- clinical_filter$EXAMDATE
clinical_filter$EXAMDATE2[!is.na(diagnosis_filter$DX2)] <- diagnosis_filter$EXAMDATE[!is.na(diagnosis_filter$DX2)]

# mutate DX with non-Dementia & Dementia samples at last follow up
pheno_baseline_cn$DX2 <- clinical_filter$DX2[match(pheno_baseline_cn$RID, clinical_filter$RID)]
pheno_baseline_cn$EXAMDATE2 <- clinical_filter$EXAMDATE2[match(pheno_baseline_cn$RID, clinical_filter$RID)]
survdate <- interval (
  pheno_baseline_cn$EXAMDATE,
  pheno_baseline_cn$EXAMDATE2) %>% 
  time_length(unit = "days")
pheno_baseline_cn$SurvDays <- survdate

# match to adni_se
adni_se_surv <- adni.se[,match(pheno_baseline_cn$RID, adni.se$RID)]
colData(adni_se_surv) <- DataFrame(pheno_baseline_cn)
colnames(adni_se_surv) <- adni_se_surv$barcodes

saveRDS(
  adni_se_surv,
  file.path(dir.data.pca, "ADNI_all_sample_cn_mci_baseline_surv.rds")
)

