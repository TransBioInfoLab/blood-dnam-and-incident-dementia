#######################################################################################################
# =================================================================================================== #
# Function for annotation and inflation adjusted 
# =================================================================================================== #
#######################################################################################################
library(SummarizedExperiment)
library(tidyverse)
library(bacon)
library(GWASTools)
library(minfi)
library(rGREAT)
# ===================================================================================================
# Annotation
# ===================================================================================================
annotate_results <- function(result, 
                             array = "HM450", 
                             dir.data.aux = dir.data.aux, 
                             save = T,
                             dir.save = NULL,
                             prefix = "Framingham"){
  
  # load(file.path(dir.data.aux,"E073_15_coreMarks_segments.rda"))
  # load(file.path(dir.data.aux,"meta_analysis_cpgs.rda"))

  if(array == "HM450"){
    anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } 
  if(array == "EPIC"){
    anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }

  anno.gr <- anno %>% makeGRangesFromDataFrame(
    start.field = "pos", end.field = "pos", keep.extra.columns = T
  )
  
  # Add annotation
  result <- cbind(
    result,
    as.data.frame(anno.gr[result$cpg])[,c(1:4)] 
  )

  if(array == "HM450"){
    load(file.path(dir.data.aux,"great_HM450_array_annotation.rda"))
    result$Islands.UCSC.Relation_to_Island <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC[result$cpg,"Relation_to_Island"]
    result$UCSC_RefGene_Name <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other[result$cpg,"UCSC_RefGene_Name"]       
    result$UCSC_RefGene_Group <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other[result$cpg,"UCSC_RefGene_Group"]     
  }
  if(array == "EPIC"){
    load(file.path(dir.data.aux,"great_EPIC_array_annotation.rda"))
    result$Islands.UCSC.Relation_to_Island <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC[result$cpg,"Relation_to_Island"]
    result$UCSC_RefGene_Name <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[result$cpg,"UCSC_RefGene_Name"]       
    result$UCSC_RefGene_Group <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[result$cpg,"UCSC_RefGene_Group"]  
  }
  result <- dplyr::left_join(result, great, by = c("seqnames","start","end","cpg"))
  
  if(save){
    write_csv(
      result,
      file.path(dir.save, paste0(prefix, "_annotated_results.csv"))
    )
  }
  
  return(result)
}

add_chmm_annotation <- function(result, dir.data.aux = dir.data.aux) {
  
  load(file.path(dir.data.aux,"E073_15_coreMarks_segments.rda"))
  
  message("Annotating E073_15_coreMarks_segments")

  result$region <- paste0(result$seqnames,":",result$start,"-", result$end)      
  result$start <- as.numeric(result$start)
  result$end <- as.numeric(result$end)
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )
  hits <- findOverlaps(result.gr, ChmmModels.gr) %>% as.data.frame()
  hits$state <- ChmmModels.gr$state[hits$subjectHits]
  hits$region <- result$region[hits$queryHits]
  result$E073_15_coreMarks_segments_state <- hits$state[match(result$region,hits$region)]
  result$region <- NULL
  
  result
}

# Combp annotation
add_combp_annotation <- function(result, cpg, array = "EPIC", dir.data.aux){
  
  # load auxilary data
  load(file.path(dir.data.aux,"E073_15_coreMarks_segments.rda"))
  
  data <- readr::read_tsv(
    file.path(dir.data.aux,"AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
  )
  CellType.selected <- readxl::read_xlsx(
    file.path(dir.data.aux,"Nassser study selected biosamples.xlsx"),col_names = FALSE
  ) %>% dplyr::pull(1)
  
  data.filtered <- data %>% dplyr::filter(CellType %in% CellType.selected) %>% 
    dplyr::filter(!isSelfPromoter)  %>% 
    dplyr::filter(class != "promoter")
  
  nasser.enhancer.gr <- data.filtered %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "chr",
    keep.extra.columns = TRUE
  )
  
  message("Annotating E073_15_coreMarks_segments")
  result$seqnames <- paste0("chr",result$`#chrom`) 
  result$region <- paste0(result$seqnames,":",result$start,"-", result$end)      
  result$start <- as.numeric(result$start)
  result$end <- as.numeric(result$end)
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )
  hits <- findOverlaps(result.gr, ChmmModels.gr) %>% as.data.frame()
  hits$state <- ChmmModels.gr$state[hits$subjectHits]
  hits$region <- result$region[hits$queryHits]
  result$E073_15_coreMarks_segments_state <- hits$state[match(result$region,hits$region)]
  
  message("Annotating GREAT")
  
  job <- submitGreatJob(result.gr, species = "hg19")
  regionsToGenes_gr <- rGREAT::getRegionGeneAssociations(job)
  regionsToGenes <- as.data.frame(regionsToGenes_gr)

  GREAT_annotation <- lapply(
    seq_len(length(regionsToGenes$annotated_genes)),
    function(i) {
      g <- ifelse(regionsToGenes$dist_to_TSS[[i]] > 0,
                  paste0(regionsToGenes$annotated_genes[[i]],
                         " (+", regionsToGenes$dist_to_TSS[[i]], ")"),
                  paste0(regionsToGenes$annotated_genes[[i]],
                         " (", regionsToGenes$dist_to_TSS[[i]], ")"))
      paste0(g, collapse = ";")
    }
  )
  great <-  dplyr::select(
    regionsToGenes, "seqnames", "start", "end", "width")
  great <- data.frame(
    great, GREAT_annotation = unlist(GREAT_annotation))
  
  
  result <- dplyr::left_join(result, great,by = c("seqnames","start","end"))
  
  message("Annotating Island")
  
  result$chrom <- result$seqnames
  result <- coMethDMR::AnnotateResults(result, arrayType = array, nCores_int = 10)
  result$chrom <-  result$chr <- NULL
  
  message("Annotating enhancer")
  
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )
  hits <- findOverlaps(result.gr, nasser.enhancer.gr) %>% as.data.frame()
  result$nasser_is_enahncer <- FALSE
  result$nasser_is_enahncer[unique(hits$queryHits)] <- TRUE
  result$nasser_is_enahncer_cell_types <- NA
  result$nasser_is_enahncer_cell_types[unique(hits$queryHits)]  <- sapply(unique(hits$queryHits), function(x){nasser.enhancer.gr$CellType[hits$subjectHits[hits$queryHits %in% x]] %>% unique %>% paste(collapse = ",")})
  
  result$cpgs_in_region <- sapply(
    result$region, 
    function(x){
      coMethDMR:::GetCpGsInRegion(x,genome = "hg19",arrayType = array) %>% 
        intersect(cpg)%>% 
        paste(collapse = ",")
    })
  return(result)
  
}
# ===================================================================================================
# Bacon inflation
# ===================================================================================================
bacon_adj <- function(data, est_var, z_var, std_var, 
                      use_z = F, 
                      save = F, 
                      dir.save = NULL,
                      prefix = "Framingham"){
  
  ### 1. Compute genomic inflation factor before bacon adjustment
  data <- data %>% mutate(
    chisq = get(z_var)^2
  ) 
  
  # inflation factor - last term is median from chisq distrn with 1 df  
  inflationFactor <- median(data$chisq,na.rm = TRUE) / qchisq(0.5, 1)
  print("lambda")
  print(inflationFactor)

  ### 2. bacon analysis
  if(use_z){
    ### 2. bacon analysis
    z_scores <- data[[z_var]]
    bc <- bacon(
      teststatistics = z_scores,
      na.exclude = TRUE,
      verbose = F
    )
    # inflation factor
    print("lambda.bacon")
    print(inflation(bc))
    # bias
    print("estimate bias")
    print(bias(bc))
    print("estimates")
    print(estimates(bc))
    
    ### 3. Create final dataset
    data.with.inflation <- data %>% mutate(
      zScore.bacon = tstat(bc)[,1],
      pValue.bacon.z = pval(bc)[,1],
      fdr.bacon.z = p.adjust(pval(bc), method = "fdr"),
    ) %>% mutate(z.value = z_scores)
    
    print("o After bacon correction")
    print("Conventional lambda")
    lambda.con <- median((data.with.inflation$zScore.bacon) ^ 2,na.rm = TRUE)/qchisq(0.5, 1)
    print(lambda.con)
    
    percent_null <- trunc ( estimates(bc)[1]*100, digits = 0)
    percent_1  <- trunc ( estimates(bc)[2]*100, digits = 0 )
    percent_2  <- 100 - percent_null - percent_1  
    bc2 <- bacon(
      teststatistics = data.with.inflation$zScore.bacon,
      na.exclude = TRUE,
      priors = list(
        sigma = list(alpha = 1.28,  beta = 0.36), 
        mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
        epsilon = list(gamma = c(90, 5, 5)))
    )
  } else {
    est <- data[[est_var]]
    se <- data[[std_var]]
    
    
    bc <- bacon(
      teststatistics = NULL,
      effectsizes = est,
      standarderrors = se,
      na.exclude = TRUE,
      verbose = F
    )
    #  posteriors(bc)
    # inflation factor
    print("lambda.bacon")
    print(inflation(bc))
    # bias
    print("estimate bias")
    print(bias(bc))
    print("estimates")
    print(estimates(bc))
    
    ### 3. Create final dataset
    data.with.inflation <- data.frame(
      data,
      Estimate.bacon = bacon::es(bc),
      StdErr.bacon = bacon::se(bc),
      pValue.bacon = pval(bc),
      fdr.bacon = p.adjust(pval(bc), method = "fdr"),
      stringsAsFactors = FALSE
    )
    
    print("o After bacon correction")
    print("Conventional lambda")
    lambda.con <- median((data.with.inflation$Estimate.bacon/data.with.inflation$StdErr.bacon) ^ 2,na.rm = TRUE)/qchisq(0.5, 1)
    print(lambda.con)
    
    # percent_null <- trunc ( estimates(bc)[1]*100, digits = 0)
    # percent_1  <- trunc ( estimates(bc)[2]*100, digits = 0 )
    # percent_2  <- 100 - percent_null - percent_1  
    bc2 <- bacon(
      teststatistics = NULL,
      effectsizes =  data.with.inflation$Estimate.bacon,
      standarderrors = data.with.inflation$StdErr.bacon,
      na.exclude = TRUE,
      priors = list(
        sigma = list(alpha = 1.28,  beta = 0.36), 
        mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
        epsilon = list(gamma = c(99, .5, .5)))
    )
  }
  
  print("inflation")
  print(inflation(bc2))
  print("estimates")
  print(estimates(bc2))
  data.with.inflation$chisq <- NULL
  
  inflation.stat <- data.frame(
    "Inflation.org" = inflationFactor,
    "Inflation.bacon" = inflation(bc),
    "Bias.bacon" = bias(bc),
    "Inflation.after.correction" = lambda.con,
    "Inflation.bacon.after.correction" = inflation(bc2),
    "Bias.bacon.after.correction" = bias(bc2)
  )
  
  
  if(save){
    readr::write_csv(
      data.with.inflation,
      file.path(dir.save, paste0(prefix, "_bacon_correction.csv"))
    )
    
    writexl::write_xlsx(
      inflation.stat,
      file.path(dir.save, paste0(prefix, "_inflation_stats.xlsx"))
    )
  }
  return(
    list(
      "data.with.inflation" = data.with.inflation,
      "bacon.obj" = bc,
      "inflation.stat" = inflation.stat
    )
  )
}



