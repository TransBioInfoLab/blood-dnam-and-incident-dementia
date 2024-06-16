dir.base <- "."
dir.data <- file.path(dir.base,"../DATASETS/", cohort)
dir.data.combp <- file.path(dir.base, "datasets/combp")
dir.data.processed <- file.path(dir.data, "DNAm/processed/EPIC") 
dir.results <- file.path(dir.base, "analysis_results")
dir.results.combp <- file.path(dir.results, "combp")
dir.results.meta <- file.path(dir.results, "meta_analysis")
dir.data.aux <- file.path(dir.base, "../DATASETS/Aux_Sync")
for(p in grep("dir.",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

# Auxilary function
source(file.path(dir.base, "code_validation/utility/annotation_and_bacon.R"))

# Load data

results <- read_csv(
  file.path(dir.results.meta,
            "meta_analysis_FHS9_ADNI_results.csv")
)

combp_results <- readxl::read_xlsx(
  file.path(dir.data.combp, "cnew.regions-p.bed_4-19-2024.xlsx")
)

# Add annotation
results_anno <- add_combp_annotation(
  result = combp_results,
  cpg = results$cpg,
  array = "EPIC",
  dir.data.aux = dir.data.aux
)

# Add direction
direction <- plyr::laply(
  str_split(results_anno$cpgs_in_region, ","),
  .fun = function(cpgs){
    est <- results %>% filter(cpg %in% cpgs) %>% pull(estimate)
    paste0(ifelse(est < 0, "-", "+"), collapse = "")
  }
)
results_anno$direction <- direction
# Calculate percentage of direction
pct.direction <- plyr::laply(
  str_split(direction, ""),
  .fun = function(d){
    sum(d == "-")/length(d)
  }
)
results_anno$pct_direction <- pct.direction

write_csv(
  results_anno,
  file.path(dir.results.combp, "combp_FHS9_ADNI_results_annotated.csv")
)

# Check with pedno
framingham_se <- readRDS(
  file.path(dir.data.processed, "FHS_OffSpring_EPIC_rm_dem_nw_se.RDS")
)
## Load meta results
library(lme4)
cpgs.sig <- read_csv(
  file.path(dir.results.meta,
            "meta_analysis_FHS9_ADNI_results.csv")
) %>% filter(pVal.fixed < 1e-05 & 
               direction %in% c("++", "--") & 
               ADNI_pValue.bacon < 0.05 & FHS9_pValue.bacon < 0.05) %>%
  pull(cpg) ## 45

## Load comb-p DMR results
sig.dmr <- read_csv(
  file.path(dir.results.combp, "combp_FHS9_ADNI_results_annotated.csv")
) %>% 
  filter(z_sidak_p < 0.05 & z_p < 1E-05 & pct_direction %in% c(0,1) & n_probes >= 3)
sig.dmr <- str_split(sig.dmr$cpgs_in_region, ",") %>% unlist() ## 233

betaMat <- assay(framingham_se)
pheno <- colData(framingham_se) %>% data.frame()
M <- minfi::logit2(betaMat)
M <- M[unique(c(cpgs.sig, sig.dmr)),]

results_mm <- plyr::adply(
  M,
  1,
  .fun = function(mvalue){
    mod = lmer(mvalue ~ (1|pedno), data.frame(mvalue = mvalue, pheno))
    vc <- VarCorr(mod) %>% as.data.frame()
    data.frame(ICC = vc$vcov[1]/(vc$vcov[1] + vc$vcov[2]))
  },.id = "cpg"
)

write_csv(
  results_mm,
  file.path(dir.results.meta, "ICC_pedno_sig_cpgs_meta_dmr.csv")
)