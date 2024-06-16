#######################################################################################################
# =================================================================================================== #
# Preprocessing
# =================================================================================================== #
#######################################################################################################
# Load packages 
library(data.table)
library(tidyverse)
library(SummarizedExperiment)
library(ExperimentHub)
library(minfi)
library(wateRmelon)
library(DMRcate)
library(readxl)
library(RPMM)
library(sm)
library(S4Vectors)
library(survival)
library(ggplot2)
library(ggrepel)
library(EpiDISH)
# ---------------------------------------------------------------------------------------------------
# Create dir 
# ---------------------------------------------------------------------------------------------------
create_dir <- function(cohort, dir.base, create_pheno_dir = T){

  dir.data <<- file.path(dir.base,"../DATASETS/", cohort)
  dir.data.raw <<- file.path(dir.data, "DNAm/raw/") 
  dir.data.processed <<- file.path(dir.data, "DNAm/processed/") 
  dir.data.pca <<- file.path(dir.data, "DNAm/pca_filtering/") 
  if(create_pheno_dir){
    dir.data.pheno.raw <<- file.path(dir.data, "Phenotype/raw/")
    dir.data.clinical <<- file.path(dir.data, "Phenotype/processed/")
  }
  # Create results dir
  dir.results <<- file.path(dir.base, "analysis_results/", cohort)
  
}
#####################################################################################################
# Sample QC
#####################################################################################################
# ---------------------------------------------------------------------------------------------------
# Filter bs  
# ---------------------------------------------------------------------------------------------------
samples_filter_bs <- function(RGSet, threshold = 85, save = F, dir_save = NULL){
  
  # Calculate sample bisulfiteConversion
  bs <- data.frame(bisulfiteConversion = bscon(RGSet))
  
  # Filter bs > threshold
  bsFilteredOut <- row.names(bs)[bs$bisulfiteConversion < threshold]
  RGSet_filtered_bs <- RGSet[,!colnames(RGSet) %in% bsFilteredOut]
  
  # Return a list of filtered RGChannelSet and bs
  out <- list(
    RGSet = RGSet_filtered_bs,
    bs = bs
  )
  
  nb_samples_after_bs_filtered <<- ncol(RGSet_filtered_bs)
  
  if(save){
    # If save the results, save an rda file to the predefined saving path dir_save
    save(
      RGSet_filtered_bs,
      bs,
      file = file.path(dir_save, paste0("/RGSet_filtered_bs_min_", threshold,".rda"))
    )
  }
  
  return(out)
}
# ---------------------------------------------------------------------------------------------------
# Filter Unmatch sex
# ---------------------------------------------------------------------------------------------------
samples_filter_sex <- function(object, sex_col, clinical_df = NULL, 
                              package = "minfi", check_samples = T, verbose = F){
  
  # Two package can be used to predict sex: minfi and wateRmelon
  # Default is minfi
  # The Sex column should be male and female, both upper and lower letters are good
  
  if(package == "minfi"){
    if(!class(object)[1] %in% c("MethylSet", "RGChannelSet")){
      stop("Please convert object to MethylSet or RGChannelSet")
    } else {
      # Convert to GenomicRatioSet
      Gset <- mapToGenome(object)
      
      # Predicted gender with minfi
      pred_sex <- minfi::getSex(Gset)
      pred_sex <- ifelse(pred_sex$predictedSex == "M", "MALE", "FEMALE")
    }
    
    if(is.null(clinical_df)){
      clinical_df <- colData(object)
    }
  }
  
  if(package == "wateRmelon"){
    if(class(object)[1] %in% c("MethylSet", "RGChannelSet")){
      beta <- getBeta(object)
    } 
    
    pred_XY <- wateRmelon::estimateSex(beta, do_plot = F)
    pred_sex <- ifelse(pred_XY$predicted_sex == "Male", "MALE", "FEMALE")
  }
  
  record_sex <- toupper(clinical_df[[sex_col]])
  if(verbose){
    cat("Identical with recorded sex: ", identical(pred_sex, record_sex), "\n")
    cat("Index of not identical samples:", which(pred_sex != record_sex))
  }
  
  # Remove samples with unmatched sex
  idx <- which(pred_sex == record_sex)
  object <- object[,idx]
  clinical_df <- clinical_df[idx,]
  
  if(check_samples){
    # Assign number counting into a vector 
    nb_samples_um_sex <<- ncol(object)
  }
  
  if(is.null(clinical_df)){
    out <- object
  } else {
    out <- list(
      object = object,
      clinical_df = clinical_df,
      rm_idx = which(pred_sex != record_sex)
    )
  }
  
  return(out)
}
# ---------------------------------------------------------------------------------------------------
# Wrapper function
# ---------------------------------------------------------------------------------------------------
samples_filter_qc <- function(RGSet, threshold = 85, sex_col, 
                              clinical_df = NULL, package = "minfi", 
                              check_samples = T, save = F, dir_save = NULL, verbose = F){
  
  # Filter bs  
  RGSet <- samples_filter_bs(
    RGSet = RGSet, threshold = threshold, save = save, dir_save = dir_save
  )
  
  if(!is.null(clinical_df)){
    clinical_df <- clinical_df[colnames(RGSet$RGSet),]
  }
  
  # Filter Unmatch sex
  RGSet_filtered <- samples_filter_sex(
    RGSet = RGSet$RGSet,
    sex_col = sex_col,
    clinical_df = clinical_df,
    package = package, 
    check_samples = check_samples,
    verbose = verbose
  )
  
  return(RGSet_filtered)
}
# ---------------------------------------------------------------------------------------------------
# Filter outliers by PCA
# ---------------------------------------------------------------------------------------------------
samples_filter_pca <- function(object, save = F, dir_save = NULL, 
                               check_samples = T, return_model = T,...){
  
  # Object be a SummarizedExperiment set
  exp_mat <- assay(object)
  
  # order by most variable probes on top
  betaOrd_mat <- OrderDataBySd(exp_mat)
  
  # Select top 50000 variable probes and run PCA
  pca <- prcomp(
    t(betaOrd_mat[1:50000,]),
    center = TRUE,
    scale = TRUE,
    ...
  )
  
  # Select first two PCs
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  
  # Calculate 3sd outlier bands of first two PCs
  out3sd_df <- d %>% 
    summarise(across(everything(), list(lower = ~ mean(.x) - 3 * sd(.x),
                                        upper = ~ mean(.x) + 3 * sd(.x))))
  
  # Find outliers by PC larger or smaller than first two PCs' bands
  d <- d %>% 
    mutate(outlier_PC1 = ifelse(PC1 >= out3sd_df$PC1_lower & PC1 <= out3sd_df$PC1_upper, 0, 1),
           outlier_PC2 = ifelse(PC2 >= out3sd_df$PC2_lower & PC2 <= out3sd_df$PC2_upper, 0, 1))
  
  # Filter outlier
  noOutliers <- d[which(d$outlier_PC1 == 0 & d$outlier_PC2 == 0), ]
  object <- object[, rownames(noOutliers)]
  
  if(save){
    # Save PCA results and outliers
    write.csv(d, file.path(dir_save, "PCs_usingBetas.csv"))
    saveRDS(pca, file.path(dir_save, "PCA_model_usingBetas.RDS"))
    saveRDS(object, 
            file.path(dir_save, "SummarizedExperiment_PCfiltered.RDS"))
  }
  
  if(check_samples){
    nb_samples_after_pca <<- ncol(object)
  }
  
  if(return_model){
    out <- list(
      SE = object,
      pca = pca,
      PC_info = d
    )
  } else {
    out <- SE
  }
  
  return(out)
  
}
# ---------------------------------------------------------------------------------------------------
## Inner function for pc
# ---------------------------------------------------------------------------------------------------
OrderDataBySd <- function(exp_mat){
  # compute sds for each row
  sds <- matrixStats::rowSds(exp_mat)
  sdsSorted <- order(sds, decreasing = TRUE)
  
  # order by most variable probes on top
  exp_mat[sdsSorted ,]
}
plotPCA <- function (pca, dataset, pheno, group_char, ntop){
  
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  # merge pheno info with PCs
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2])
  d <- d[as.character(pheno$sample) ,]
  
  test <- identical (as.character(pheno$sample), row.names (d))
  
  if (test == TRUE){
    
    plotData <- merge (d, pheno, by.x = "row.names", by.y = "sample")
    
    
    # add lines for 3 SD from mean
    meanPC1 <- mean (plotData$PC1)
    sdPC1   <- sd (plotData$PC1)
    
    meanPC2 <- mean (plotData$PC2)
    sdPC2   <- sd (plotData$PC2)
    
    # add flag for outlier samples
    plotData$outlier <- ifelse ( abs(plotData$PC1) > meanPC1 + 3*sdPC1 | abs(plotData$PC2) > meanPC2 + 3*sdPC2,
                                 1, 0 )
    
    plotData$outlier_name <- ifelse ( abs(plotData$PC1) > meanPC1 + 3*sdPC1 | abs(plotData$PC2) > meanPC2 + 3*sdPC2,
                                      plotData$Row.names, "" )
    
    title <- paste0("dataset = ", dataset, ", top = ", ntop, " probes ")
    subtitle <- paste0(" x: mean +/- 3*sdPC1 = ", round(meanPC1,1), " +/- 3*", round(sdPC1,1) ,
                       "     y: mean +/- 3*sdPC2 = ", round(meanPC2,1), " +/- 3*", round(sdPC2,1))
    plotData[[group_char]] <- as.factor(plotData[[group_char]])
    
    p <- ggplot(data= plotData, aes_string(x="PC1", y="PC2", color = group_char)) +
      geom_point(size=1) +
      theme_bw() +
      xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
      ggtitle(title, subtitle = subtitle) +
      geom_hline (yintercept = meanPC2 + 3*sdPC2, linetype = "dashed") +
      geom_hline (yintercept = meanPC2 - 3*sdPC2, linetype = "dashed") +
      
      geom_vline (xintercept = meanPC1 + 3*sdPC1, linetype = "dashed") +
      geom_vline (xintercept = meanPC1 - 3*sdPC1, linetype = "dashed") +
      geom_text_repel (aes(label = outlier_name), show.legend = FALSE)
    
    print (p)
    
    
    return (plotData)
    
  }
}
#####################################################################################################
# Converting Functions
#####################################################################################################
# ---------------------------------------------------------------------------------------------------
# Convert to SummarizedExperiment
# ---------------------------------------------------------------------------------------------------
convert_SE <- function(object, pheno_df = NULL, array = "HM450", genome = "hg19"){
  
  # Check the object's class 
  if(class(object)[1] %in% c("MethylSet", "RGChannelSet")){
    data <- getBeta(object)
  } else {data <- object}
  
  # Check if the pheno data is provided
  if(!is.null(pheno_df)){
    coldata <- DataFrame(pheno_df)
  } else {
    coldata = colData(object)
  }
  
  # Check the array type
  if(array == "HM450"){
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else {
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }

  # Get annotation
  anno.gr <- anno %>% makeGRangesFromDataFrame(
    start.field = "pos", end.field = "pos", keep.extra.columns = T
  )
  
  # Add annotation
  data <- data[rownames(data) %in% names(anno.gr),]
  rowData <- anno.gr[rownames(data),]
  data <- as.matrix(data)
  
  # Convert to SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = SimpleList("DNAm" = data),
    rowData = rowData,
    metadata = list("genome" = genome),
    colData = coldata
  )
  
  return(se)
}
# ---------------------------------------------------------------------------------------------------
# Convert to MethylSet
# ---------------------------------------------------------------------------------------------------
convert_MethylSet <- function(RGSet, probes_keep = NULL, save = F, dir_save = NULL){
  
  # Create MethylSet from RGChannelSet
  MSet <- preprocessRaw(RGSet)
  
  if(!is.null(probes_keep)){
    # Filter MethylSet by probes_keep
    MSet <- MSet[rownames(MSet) %in% probes_keep,]
  }
  
  if(save){
    save(
      MSet, 
      file =  file.path(dir_save,"MethylSet.rda")
    )
  }
  
  return(MSet)
  
}
#####################################################################################################
# Probes QC
#####################################################################################################
# ---------------------------------------------------------------------------------------------------
# Filter detection p
# ---------------------------------------------------------------------------------------------------
probes_filter_p <- function(RGSet, p_thres = 0.05, sample_prop = 0.9, 
                            return_beta = T, save = F, dir_save = NULL){
  
  # Calculate detection P
  detP <- detectionP(RGSet, type = "m+u")
  
  # Filter by P < p_thres in >= sample_prop
  failed.01 <- detP > p_thres
  sample_thres <- 1 - sample_prop
  if(sample_thres == 0){
    passedProbes <- rownames(failed.01)[rowMeans(failed.01) == 0] 
  } else {
    passedProbes <- rownames(failed.01)[rowMeans(failed.01) < sample_thres] 
  }
  
  if(save){
    # If save is true, save the passedProbes for future use
    save(passedProbes,
         file = file.path(dir_save,"detectionP_passed_probes.rda"))
  }
  
  if(return_beta){
    # Get beta matrix from RGSet
    betaSet <- getBeta(RGSet)
    
    # Select probes pass detection P
    beta_filtered <- betaSet[rownames(betaSet) %in% passedProbes, ]
    
    out <- list(
      beta = beta_filtered,
      passedProbes = passedProbes
    )
  } else {
    out <- passedProbes
  }
  
  return(out)
  
}
# ---------------------------------------------------------------------------------------------------
# Filter CG, SNP with MAP and XY
# ---------------------------------------------------------------------------------------------------
probes_filter_SNP <- function(object, rmNonCG = T, 
                              rmSNP = T,
                              rmcrosshyb = T, rmXY = T, 
                              dist = 5, mafcut = 0.01, and = T, check_probes = T, 
                              save = F, dir_save = NULL){
  # Object should be a beta or M matrix
  
  # Remove non-cg probes
  if(rmNonCG){
    object <- object[grep("cg",rownames(object)),]
  }

  # Filter probes by SNP with MAF > mafcut was present in the last dist bp of the probe and XY probes
  if(rmSNP) {
    beta_filtered <- rmSNPandCH(
      object = object,
      dist = dist,
      mafcut = mafcut,
      and = and,
      rmcrosshyb = rmcrosshyb,
      rmXY = rmXY
    )
  } else {
    beta_filtered <- object
  }
  
  if(check_probes){
    num_probes <- c(
      "nb_probes_after_remove_noncg" = nrow(object),
      "nb_probes_after_max_dist_mafcut" = nrow(beta_filtered)
    )
    
    out <- list(
      beta = beta_filtered,
      num_probes = num_probes
    )
  } else {
    out <- beta_filtered
  }
  
  if(save){
    saveRDS(
      out,
      file = file.path(dir_save,"BetaSet_after_QC.rds")
    )
  }
  
  return(out)
}
# ---------------------------------------------------------------------------------------------------
# Wrapper function
# ---------------------------------------------------------------------------------------------------
probes_filter_qc <- function(RGSet, p_thres = 0.05, sample_prop = 0.9, 
                             rmNonCG = T, rmSNP = T, rmXY = T, dist = 5, 
                             mafcut = 0.01, and = T, convert_MethylSet = T,
                             check_probes = T, save = F, dir_save = NULL){
    
    # Filter detection p and get beta
    p_list <- probes_filter_p(RGSet, 
                              p_thres = p_thres, 
                              sample_prop = sample_prop, 
                              save = save, 
                              dir_save = dir_save)
    beta <- p_list$beta
    
    # Filter cg, SNP with MAP and XY
    beta_ls <- probes_filter_SNP(
      object = beta,
      rmNonCG = rmNonCG,
      rmSNP = rmSNP,
      rmXY = rmXY,
      dist = dist,
      mafcut = mafcut, 
      and = and,
      check_probes = check_probes,
      save = save,
      dir_save = dir_save,
    )
    
    # Check probes
    if(check_probes){
      nb_probes_count <<- c("nb_probes_after_p" = nrow(beta), 
                           beta_ls$num_probes)
      betaSet <- beta_ls$beta
    } else {betaSet <- beta_ls}
    
    # Convert MethylSet
    if(convert_MethylSet){
      MSet <- convert_MethylSet(
        RGSet = RGSet,
        probes_keep = rownames(betaSet),
        save = save,
        dir_save = dir_save
      )
      return(MSet)
    } else {
      return(betaSet)
    }
  
  }
#####################################################################################################
# Normalization
#####################################################################################################
methyl_norm <- function(object, method = "dasen", save = F, 
                        dir_save = NULL, array = "HM450", parallel = F){
  
  # Two method can be selected: dasen and BMIQ
  
  if(method == "dasen"){
    if(!class(object)[1] %in% c("MethylSet")){
      stop("Please convert object to MethylSet")
    } else {
      norm_set <- dasen(object)
      
      if(save){
        saveRDS(norm_set, file = file.path(dir_save,"MethylSet_dasen.rds"))
      }
    }
  }
  
  if(method == "BMIQ"){
    betaQN <- lumiN(x.lumi = object, method = "quantile")
    
    if(array == "HM450"){
      annotType <- 
        IlluminaHumanMethylation450kanno.ilmn12.hg19::Manifest
    } else {
      annotType <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest
    }
    annotType$designTypeNumeric <- ifelse(annotType$Type == "I",1,2)
    type12 <- annotType$designTypeNumeric[match(rownames(betaQN),rownames(annotType))]
    
    set.seed (946)
    if(parallel){
      doParallel::registerDoParallel(cores = 4)
    }
    norm_set <- plyr::aaply(
      betaQN, 2,
      function(x){
        norm_ls <- BMIQ(x, design.v = type12, plots = FALSE)
        return (norm_ls$nbeta)
      },.progress = "time",.parallel = parallel
    ) %>% t()
    
    if(save){
      saveRDS(norm_set, file.path(dir_save,"BetaSet_QNBMIQ.RDS"))
    }
  }
  
  return(norm_set)
}

#####################################################################################################
# The whole process wrapper function
#####################################################################################################
methyl_qc <- function(RGSet, clinical = NULL, 
                     # Filter bs  
                     bs_thres = 85,
                     # detectionP
                     p_thres = 0.05, sample_prop = 0.9, 
                     # Filter CG, SNP with MAP and XY
                     rmNonCG = T, rmSNP = T, rmXY = T, 
                     dist = 5, mafcut = 0.01, and = T,
                     # Filter Unmatch sex
                     pred_sex_package = "minfi", sex_col,
                     # Normalization
                     method_norm = "dasen", array = "HM450", parallel = F,
                     # Outlier detection
                     genome = "hg19",
                     # Saving options
                     check_steps = T, save = T, create_pheno_dir = F,
                     dir = getwd(), cohort_name = "FHS"){
  
  # Create dir
  if(save){
    if(!is.null(clinical)) create_pheno_dir <- T
    create_dir(cohort = cohort_name, 
               dir.base = dir, 
               create_pheno_dir = create_pheno_dir)
    for(p in grep("dir.",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)
  }
  
  nb_probes_before_qc <<- nrow(RGSet)
  nb_samples_before_qc <<- ncol(RGSet)
  
  # Sample QC
  RGSet_filtered <- samples_filter_qc(RGSet, 
                                      threshold = bs_thres, 
                                      sex_col = sex_col, 
                                      clinical_df = clinical, 
                                      package = package, 
                                      check_samples = check_steps,
                                      save = save, 
                                      dir_save = dir.data.processed)
  if(!is.null(clinical)){
    clinical <- RGSet_filtered$clinical_df
    RGSet_filtered <- RGSet_filtered$object
  }
  
  # Probes QC,
  Mset <- probes_filter_qc(RGSet_filtered,
                           p_thres = p_thres, sample_prop = sample_prop, 
                           rmNonCG = rmNonCG, rmSNP = rmSNP, rmXY = rmXY, dist = dist, 
                           mafcut = mafcut, and = and, 
                           check_probes = check_steps,
                           save = save, dir_save = dir.data.processed,
                           convert_MethylSet = T)
  
  # Normalization
  if(method_norm == "BMIQ"){
    object <- getBeta(Mset)
  }
  if(method_norm == "dasen"){
    object = Mset
  }
  norm_set <- methyl_norm(
    object = object, 
    method = method_norm,
    save = save, 
    dir_save = dir.data.processed,
    parallel = parallel,
    array = array
  )

  # Convert to SummarizedExperiment
  SE <- convert_SE(
    object = norm_set,
    pheno_df = clinical,
    array = array, 
    genome = genome
  )
  
  # Filter outliers by PCA
  SE_final <- samples_filter_pca(
    object = SE, save = save, dir_save = dir.data.pca, 
    check_samples = check_steps, return_model = F,
  )
  
  if(check_steps){
    # Samples check
    df_samples <- data.frame(Description = c(
      "Number of Samples before QC", 
      "Number of Samples after filtered bs" ,
      "Number of Samples after filtered unmatched sex",
      "Number of Samples after filtered outliers by PCA"
    ),
    Count = c(
      nb_samples_before_qc, 
      nb_samples_after_bs_filtered,
      nb_samples_um_sex,
      nb_samples_after_pca
    ))
    
    # Probes check
    df_probes <- data.frame(
      Description = c("Number of Probes before QC",
                      "Number of Probes after filtered detection P",
                      "Number of Probes after filtered non-cg",
                      "Number of Probes with no crosshyb, no X, no Y, mafcut = 0.01"), 
      Count = c(nb_probes_before_qc, nb_probes_count))
    
    step_count <<- list(sample_steps = df_samples, probe_steps = df_probes)
  }
  
  if(save){
    write_xlsx(
      step_count,
      file.path(dir.data.pca, paste0(cohort_name, "_Step_count.xlsx"))
    )
    
    if(!is.null(clinical)) {
      clinical_df <- clinical[colnames(SE_final),]
      write_csv(
        clinical_df,
        file.path(dir.data.clinical, paste0(cohort_name, "_clinical_data.csv"))
      )
    }
  }
  
  return(SE_final)
}
