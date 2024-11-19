# -----------------------------------------------------------------------------------------------------------
# For reproducible research, please install the following R packages 
# and make sure the R and BiocManager versions are correct
# Session Info ----------------------------------------------------------------------------------------------
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       macOS Sonoma 14.5
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-06-03
# rstudio  2024.04.1+748 Chocolate Cosmos (desktop)
# pandoc   3.1.11 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64/ (via rmarkdown)
# -----------------------------------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.18")

list.of.packages <- c(
  "bacon",
  "BEclear",
  "coMethDMR",
  "data.table",
  "devtools",
  "DMRcate",
  "doParallel",
  "dorothea",
  "EpiDISH",
  "ExperimentHub",                                
  "GEOquery",                                     
  "ggpubr",    
  "GWASTools",  
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "janitor",
  "lme4",
  "matrixStats",
  "meta",
  "metafor",
  "methylGSA",
  "MethReg",
  "minfi",
  "msigdbr",
  "plyr",   
  "pROC",
  "readxl",
  "rGREAT",
  "RPMM",
  "sesame",
  "sesameData",
  "stats",                                        
  "SummarizedExperiment",    
  "survminer",
  "survival",
  "S4Vectors",
  "tidyverse",   
  "wateRmelon",
  "writexl"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  for(new in new.packages){
    if(new %in% available.packages()[,1]){
      install.packages(new)
    } else BiocManager::install(new)
  }
} 

