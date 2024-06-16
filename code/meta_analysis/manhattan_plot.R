# Manhattan plot
dir.base <- "."
dir.results <- file.path(dir.base, "analysis_results")
dir.results.meta <- file.path(dir.results, "meta_analysis")
dir.data.aux <- file.path(dir.base, "../DATASETS/Aux_Sync/") 
dir.results.combp <- file.path(dir.results, "combp")
dir.plot <- file.path(dir.results, "plots")

source(file.path(dir.base, "code/utility/plot.R")) 

results <- read_csv(
  file.path(dir.results.meta,
            "meta_analysis_FHS9_ADNI_results_annotated.csv")
)

cpgs.sig <- read_csv(
  file.path(dir.results.meta,
            "meta_analysis_FHS9_ADNI_results_annotated.csv")
) %>% filter(fdr < .05 & 
               direction %in% c("++", "--") & 
               ADNI_pValue.bacon < 0.05 & FHS9_pValue.bacon < 0.05) %>%
  pull(cpg) 

sig.dmr.df <- read_csv(
  file.path(dir.results.combp, "combp_FHS9_ADNI_results_annotated.csv")
) %>% 
  filter(z_sidak_p < 0.05 & z_p < 1E-05 & pct_direction %in% c(0,1) & n_probes >= 3) %>% 
  slice_min(z_p, n = 10)
sig.dmr <- str_split(sig.dmr.df$cpgs_in_region, ",") %>% unlist() ## 233

cpgs_df <- data.frame(
  "CpG" = results$cpg,
  "pVal.final" = results$pVal.fixed,
  "fdr" = results$fdr,
  "pos" = results$start,
  "chr" = as.numeric(gsub("chr", "", results$seqnames)),
  "GREAT" = gsub("\\(.*| ", "", results$GREAT_annotation) 
) 
cpgs_df$GREAT[cpgs_df$GREAT == "RBM14-RBM4"] <- "RBM14"

annotated_cpg <- results %>% 
  filter(cpg %in% cpgs.sig) %>%
  slice_min(fdr, n = 10) %>%
  mutate(GREAT = gsub("\\;.*| ", "", GREAT_annotation),
    TSS = as.numeric(str_extract(GREAT, "(?<=\\().*(?=\\))"))) %>%
  filter(abs(TSS) < 2000 & abs(TSS) > 100) %>%
  pull(GREAT) %>%
  gsub("\\(.*| ", "",.)
annotated_cpg[annotated_cpg == "RBM14-RBM4"]  <- "RBM14"
annotated_dmr <- results %>% 
  filter(cpg %in% sig.dmr) %>%
  mutate(GREAT = gsub("\\;.*| ", "", GREAT_annotation),
         TSS = as.numeric(str_extract(GREAT, "(?<=\\().*(?=\\))"))) %>%
  filter(abs(TSS) < 2000 & abs(TSS) > 100) %>%
  pull(GREAT) %>%
  gsub("\\(.*| ", "",.)

plot_manh(cpgs_df, annotated = unique(c(annotated_cpg, annotated_dmr)), colored = cpgs.sig)

ggplot2::ggsave(
  filename = file.path(dir.plot,"manhattan_plot.pdf"),
  width = 10,
  height = 6
)

ggplot2::ggsave(
  filename = file.path(dir.plot,"manhattan_plot.png"),
  width = 16,
  height = 9
)

