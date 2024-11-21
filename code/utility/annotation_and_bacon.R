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
# Annotate methylation results with auxiliary data
annotate_results <- function(result, 
                             array = "HM450", 
                             dir.data.aux = dir.data.aux, 
                             save = T,
                             dir.save = NULL,
                             prefix = "Framingham"){
  
  # Load the appropriate annotation data based on the array type
  if(array == "HM450"){
    anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } 
  if(array == "EPIC"){
    anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }

  # Convert annotation data into GRanges format for genomic interval operations
  anno.gr <- anno %>% makeGRangesFromDataFrame(
    start.field = "pos", end.field = "pos", keep.extra.columns = T
  )
  
  # Add annotation information to the result dataset
  result <- cbind(
    result,
    as.data.frame(anno.gr[result$cpg])[,c(1:4)]  # Add basic genomic information
  )

  # Load additional array-specific annotation data
  if(array == "HM450"){
    load(file.path(dir.data.aux, "great_HM450_array_annotation.rda"))
    result$Islands.UCSC.Relation_to_Island <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC[result$cpg, "Relation_to_Island"]
    result$UCSC_RefGene_Name <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other[result$cpg, "UCSC_RefGene_Name"]       
    result$UCSC_RefGene_Group <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other[result$cpg, "UCSC_RefGene_Group"]     
  }
  if(array == "EPIC"){
    load(file.path(dir.data.aux, "great_EPIC_array_annotation.rda"))
    result$Islands.UCSC.Relation_to_Island <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC[result$cpg, "Relation_to_Island"]
    result$UCSC_RefGene_Name <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[result$cpg, "UCSC_RefGene_Name"]       
    result$UCSC_RefGene_Group <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[result$cpg, "UCSC_RefGene_Group"]  
  }

  # Merge result with GREAT annotations
  result <- dplyr::left_join(result, great, by = c("seqnames", "start", "end", "cpg"))
  
  # Save the annotated results if specified
  if(save){
    write_csv(
      result,
      file.path(dir.save, paste0(prefix, "_annotated_results.csv"))
    )
  }
  
  return(result)
}

# Add ChromHMM segment annotations to the results
add_chmm_annotation <- function(result, dir.data.aux = dir.data.aux) {
  
  # Load the ChromHMM segments data
  load(file.path(dir.data.aux, "E073_15_coreMarks_segments.rda"))
  
  # Prepare the result as GRanges for overlap analysis
  result$region <- paste0(result$seqnames, ":", result$start, "-", result$end)      
  result$start <- as.numeric(result$start)
  result$end <- as.numeric(result$end)
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )

  # Find overlapping ChromHMM segments and annotate
  hits <- findOverlaps(result.gr, ChmmModels.gr) %>% as.data.frame()
  hits$state <- ChmmModels.gr$state[hits$subjectHits]
  hits$region <- result$region[hits$queryHits]
  result$E073_15_coreMarks_segments_state <- hits$state[match(result$region, hits$region)]
  result$region <- NULL
  
  return(result)
}

# Add enhancer and comb-p annotations to the results
add_combp_annotation <- function(result, cpg, array = "EPIC", dir.data.aux){
  
  # Load auxiliary data for enhancer annotations
  load(file.path(dir.data.aux, "E073_15_coreMarks_segments.rda"))
  data <- readr::read_tsv(
    file.path(dir.data.aux, "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
  )
  CellType.selected <- readxl::read_xlsx(
    file.path(dir.data.aux, "Nassser study selected biosamples.xlsx"), col_names = FALSE
  ) %>% dplyr::pull(1)

  # Filter enhancer data based on selected cell types
  data.filtered <- data %>% dplyr::filter(CellType %in% CellType.selected) %>% 
    dplyr::filter(!isSelfPromoter)  %>% 
    dplyr::filter(class != "promoter")
  
  nasser.enhancer.gr <- data.filtered %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "chr",
    keep.extra.columns = TRUE
  )
  
  # Annotate ChromHMM segments
  result$seqnames <- paste0("chr", result$`#chrom`) 
  result$region <- paste0(result$seqnames, ":", result$start, "-", result$end)      
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
  result$E073_15_coreMarks_segments_state <- hits$state[match(result$region, hits$region)]
  
  # Annotate with GREAT and enhancer information
  job <- submitGreatJob(result.gr, species = "hg19")
  regionsToGenes_gr <- rGREAT::getRegionGeneAssociations(job)
  regionsToGenes <- as.data.frame(regionsToGenes_gr)
  GREAT_annotation <- sapply(seq_len(nrow(regionsToGenes)), function(i) {
    dist <- regionsToGenes$dist_to_TSS[[i]]
    gene <- regionsToGenes$annotated_genes[[i]]
    paste0(gene, " (", ifelse(dist > 0, "+", ""), dist, ")", collapse = ";")
  })
  result <- dplyr::left_join(
    result, 
    data.frame(seqnames = regionsToGenes$seqnames, start = regionsToGenes$start, end = regionsToGenes$end, GREAT_annotation), 
    by = c("seqnames", "start", "end")
  )
  
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
  
  ### Step 1: Compute the genomic inflation factor before BACON adjustment
  data <- data %>% mutate(
    chisq = get(z_var)^2 # Calculate chi-square values from the z-scores
  )
  
  # Compute inflation factor (lambda) using the median chi-square statistic
  inflationFactor <- median(data$chisq, na.rm = TRUE) / qchisq(0.5, 1)
  print("lambda") # Print the inflation factor
  print(inflationFactor)

  ### Step 2: Perform BACON analysis
  if (use_z) { # Use z-scores for BACON if `use_z` is TRUE
    z_scores <- data[[z_var]]
    bc <- bacon(
      teststatistics = z_scores, # Input z-scores
      na.exclude = TRUE, # Exclude missing values
      verbose = F
    )
    
    # Print diagnostic statistics after BACON adjustment
    print("lambda.bacon")
    print(inflation(bc)) # Print the adjusted inflation factor
    print("estimate bias")
    print(bias(bc)) # Print the estimated bias
    print("estimates")
    print(estimates(bc)) # Print the estimated parameters
    
    ### Step 3: Add BACON-adjusted values to the dataset
    data.with.inflation <- data %>% mutate(
      zScore.bacon = tstat(bc)[, 1], # Adjusted z-scores
      pValue.bacon.z = pval(bc)[, 1], # Adjusted p-values
      fdr.bacon.z = p.adjust(pval(bc), method = "fdr"), # Adjusted FDR
      z.value = z_scores # Original z-scores
    )
    
    # Recalculate inflation factor after BACON adjustment
    print("o After bacon correction")
    print("Conventional lambda")
    lambda.con <- median((data.with.inflation$zScore.bacon)^2, na.rm = TRUE) / qchisq(0.5, 1)
    print(lambda.con)
    
    # Estimate null hypothesis proportions
    percent_null <- trunc(estimates(bc)[1] * 100, digits = 0)
    percent_1 <- trunc(estimates(bc)[2] * 100, digits = 0)
    percent_2 <- 100 - percent_null - percent_1
    
    # Perform a second BACON adjustment using refined priors
    bc2 <- bacon(
      teststatistics = data.with.inflation$zScore.bacon,
      na.exclude = TRUE,
      priors = list(
        sigma = list(alpha = 1.28, beta = 0.36), 
        mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
        epsilon = list(gamma = c(90, 5, 5))
      )
    )
  } else { # Use effect sizes and standard errors if `use_z` is FALSE
    est <- data[[est_var]]
    se <- data[[std_var]]
    
    bc <- bacon(
      teststatistics = NULL, # Use effect sizes and standard errors instead of z-scores
      effectsizes = est,
      standarderrors = se,
      na.exclude = TRUE,
      verbose = F
    )
    
    # Print diagnostic statistics after BACON adjustment
    print("lambda.bacon")
    print(inflation(bc))
    print("estimate bias")
    print(bias(bc))
    print("estimates")
    print(estimates(bc))
    
    ### Step 3: Add BACON-adjusted values to the dataset
    data.with.inflation <- data.frame(
      data,
      Estimate.bacon = bacon::es(bc), # Adjusted effect sizes
      StdErr.bacon = bacon::se(bc), # Adjusted standard errors
      pValue.bacon = pval(bc), # Adjusted p-values
      fdr.bacon = p.adjust(pval(bc), method = "fdr"), # Adjusted FDR
      stringsAsFactors = FALSE
    )
    
    # Recalculate inflation factor after BACON adjustment
    print("o After bacon correction")
    print("Conventional lambda")
    lambda.con <- median((data.with.inflation$Estimate.bacon / data.with.inflation$StdErr.bacon)^2, na.rm = TRUE) / qchisq(0.5, 1)
    print(lambda.con)
    
    # Perform a second BACON adjustment using refined priors
    bc2 <- bacon(
      teststatistics = NULL,
      effectsizes = data.with.inflation$Estimate.bacon,
      standarderrors = data.with.inflation$StdErr.bacon,
      na.exclude = TRUE,
      priors = list(
        sigma = list(alpha = 1.28, beta = 0.36), 
        mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
        epsilon = list(gamma = c(99, .5, .5))
      )
    )
  }
  
  # Print diagnostic statistics after BACON adjustment
  print("inflation")
  print(inflation(bc2))
  print("estimates")
  print(estimates(bc2))
  
  # Remove intermediate chi-square column from the data
  data.with.inflation$chisq <- NULL
  
  ### Step 4: Compile inflation statistics
  inflation.stat <- data.frame(
    "Inflation.org" = inflationFactor, # Original inflation factor
    "Inflation.bacon" = inflation(bc), # Adjusted inflation factor
    "Bias.bacon" = bias(bc), # Adjusted bias
    "Inflation.after.correction" = lambda.con, # Inflation factor after correction
    "Inflation.bacon.after.correction" = inflation(bc2), # After correction adjusted inflation factor
    "Bias.bacon.after.correction" = bias(bc2) # After correction adjusted bias
  )
  
  # Save results to specified directory if `save` is TRUE
  if (save) {
    readr::write_csv(
      data.with.inflation,
      file.path(dir.save, paste0(prefix, "_bacon_correction.csv"))
    )
    writexl::write_xlsx(
      inflation.stat,
      file.path(dir.save, paste0(prefix, "_inflation_stats.xlsx"))
    )
  }
  
  ### Step 5: Return adjusted data and statistics
  return(
    list(
      "data.with.inflation" = data.with.inflation,
      "bacon.obj" = bc, # BACON model object
      "inflation.stat" = inflation.stat # Inflation statistics
    )
  )
}
