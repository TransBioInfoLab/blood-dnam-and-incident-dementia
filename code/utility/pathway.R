#######################################################################################################
# ===================================================================================================
# Function for pathway analysis 
# ===================================================================================================
#######################################################################################################
library(methylGSA)
library(msigdbr)
# ===================================================================================================
# Wrapper Function for methylGSA Analysis
# ===================================================================================================
# This function performs methylGSA analysis using different methods (e.g., GSEA or ORA) on CpG data.
methylGSA_wrapper <- function(cpg.pval, array, method = "GSEA", GS, add_anno_for_sig = T, 
                              use_msigdbr = F, fun = "methylRRA", minsize = 100, maxsize = 500){
  
  # Load annotation data for significant results if needed
  if(add_anno_for_sig) {
    CpGs <- getAnnot(array.type = array) %>% as.data.frame()
  }
  
  # Loop through each gene set (GS) type and perform methylGSA analysis
  results <- plyr::llply(
    GS,
    .fun = function(gs.type) {
      
      # Optionally use msigdbr to load gene set information
      if(use_msigdbr) {
        pathway_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = toupper(gs.type))
        pathway <- unique(pathway_df$gs_name)
        
        # Create a list of gene sets
        GS.list <- plyr::llply(
          pathway,
          .fun = function(p) {
            pathway_df %>% filter(gs_name %in% p) %>% pull(gene_symbol)
          }
        )
        names(GS.list) <- pathway    
      } else {
        GS.list <- NULL
      }
      
      # Evaluate the function for methylGSA analysis
      fn <- eval(parse(text = fun))
      
      # Perform analysis with the selected function
      if(fun == "methylRRA") {
        res <- fn(
          cpg.pval = cpg.pval,
          method = method, 
          array.type = array,
          GS.type = gs.type,
          GS.list = GS.list,
          minsize = minsize,
          maxsize = maxsize
        )
      } else {
        res <- fn(
          cpg.pval = cpg.pval,
          array.type = array,
          GS.type = gs.type,
          GS.list = GS.list,
          minsize = minsize,
          maxsize = maxsize
        )
      }

      # Filter results for significant pathways (FDR < 0.05)
      res_sig <- res %>% filter(padj < 0.05)
      
      # Add annotations for significant results if specified
      if(add_anno_for_sig) {
        res_anno <- plyr::alply(
          res_sig,
          .margins = 1,
          .fun = function(pathway){

            # Determine the enrichment type (ORA or GSEA) and prepare annotations
            if(method == "ORA") {
              pathway <- pathway %>% separate_rows(overlap, sep = ",") 
              core <- "overlap"
            }
            if(method == "GSEA") {
              pathway <- pathway %>% separate_rows(core_enrichment, sep = "/") 
              core <- "core_enrichment"
            }
           
            # Annotate CpGs associated with genes in the pathway
            plyr::adply(
              pathway,
              .margins = 1,
              .fun = function(gene){
                cpgs <- CpGs %>%
                  filter(UCSC_RefGene_Name %in% gene[[core]]) %>%
                  dplyr::select(Name, UCSC_RefGene_Group)
                gene <- data.frame(gene, cpg = cpgs$Name)
              }, .id = "pathway"
            )
          }
        )
        names(res_anno) <- res_sig$Description
        
        # Return full results and annotations
        return(
          c(list(results = res),
            res_anno)
        )
      }
      
      # Return analysis results if no annotation is added
      res
    }
  )
  
  # Assign names to the results based on gene set types
  names(results) <- GS
  results
}

# ===================================================================================================
# Get Annotation Data
# ===================================================================================================
# This function retrieves annotation data for CpG probes based on the array type.
getAnnot = function(array.type, group = "all"){
  
  # Load annotation for HM450 or EPIC arrays
  if(array.type == "450K") {
    FullAnnot = tryCatch({
      getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    }, error = function(e) {
      stop("IlluminaHumanMethylation450kanno.ilmn12.hg19 needs to be installed and loaded.")
    })
  } else {
    FullAnnot = tryCatch({
      getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }, error = function(e) {
      stop("IlluminaHumanMethylationEPICanno.ilm10b4.hg19 needs to be installed and loaded.")
    })
  }
  
  # Retain relevant columns for CpG annotation
  FullAnnot = FullAnnot[, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
  
  # Filter for valid CpG IDs
  FullAnnot = FullAnnot[str_length(rownames(FullAnnot)) == 10, ]
  FullAnnot = FullAnnot[FullAnnot$UCSC_RefGene_Name != "", ]
  
  # Extract the first gene or group if multiple are listed
  FullAnnot$UCSC_RefGene_Name <- vapply(strsplit(FullAnnot$UCSC_RefGene_Name, split = ";"), '[', 1, FUN.VALUE = character(1))
  FullAnnot$UCSC_RefGene_Group <- vapply(strsplit(FullAnnot$UCSC_RefGene_Group, split = ";"), '[', 1, FUN.VALUE = character(1))
  
  # Apply filters for specific gene groups if requested
  if(group == "body") {
    FullAnnot <- FullAnnot[FullAnnot$UCSC_RefGene_Group %in% c("Body", "1stExon"), ]
  }
  
  if(group == "promoter1") {
    FullAnnot <- FullAnnot[grepl("TSS", FullAnnot$UCSC_RefGene_Group), ]
  }
  
  if(group == "promoter2") {
    FullAnnot <- FullAnnot[FullAnnot$UCSC_RefGene_Group %in% c("TSS200", "TSS1500", "1stExon", "5'UTR"), ]
  }
  
  return(FullAnnot)
}
