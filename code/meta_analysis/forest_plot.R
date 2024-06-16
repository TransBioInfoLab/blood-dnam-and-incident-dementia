# Forest plot
dir.base <- "."
dir.results <- file.path(dir.base, "analysis_results")
dir.results.meta <- file.path(dir.results, "meta_analysis")
dir.data.aux <- file.path(dir.base, "../DATASETS/Aux_Sync/") 
dir.results.combp <- file.path(dir.results, "combp")
dir.plot <- file.path(dir.results, "plots")

source(file.path(dir.base, "code/utility/plot.R")) 

# FHS exam 9
results_fhs9 <- read_csv(
  file.path(dir.results, "Framingham/EPIC_EXAM9/Framingham_EPIC_beta_single_cpg_cox_all_celltypes_bacon_results.csv")
)
results_fhs9 <- results_fhs9 %>% 
  mutate(study = "FHS",
         study_id = "FHS9")
# ADNI
results_adni <- read_csv(
  file.path(dir.results, "ADNI/cox/ADNI_adjust_celltypes_first_2pc_single_cpg_cox_bacon_correction.csv")
)
results_adni <- results_adni %>% 
  mutate(study = "ADNI",
         study_id = "ADNI")

results_meta <-  read_csv(
  file.path(dir.results.meta, "meta_analysis_FHS9_ADNI_results.csv")
)

top_cpg <- results_meta %>% 
  filter(fdr < .05,
         direction %in% c("++", "--"),
         ADNI_pValue.bacon < 0.05, 
         FHS9_pValue.bacon < 0.05) %>%
  slice_min(fdr, n = 10)

results_list <- list(
  FHS9 = results_fhs9,
  ADNI = results_adni
) %>% purrr::map(., ~.[.$cpg %in% top_cpg$cpg,])

plot_list <- plot_forest(results_list = results_list,
                         cpg_order = top_cpg$cpg,
                         effect = "Estimate.bacon",
                         se = "StdErr.bacon",
                         sm = "HR",
                         studlab = c("FHS", "ADNI"),
                         ncol = 2,
                         nrow = 5)

ggsave(
  plot = plot_list,
  height =  12,
  width = 16,
  filename = file.path(dir.plot, "Forest_plot_top_10_cpg_all_in_one.pdf")
)
