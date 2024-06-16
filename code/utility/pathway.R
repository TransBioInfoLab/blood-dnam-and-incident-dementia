#######################################################################################################
# ===================================================================================================
# Function for pathway analysis 
# ===================================================================================================
#######################################################################################################
library(methylGSA)
library(msigdbr)
# ===================================================================================================
# Wrapper function for methylGSA
# ===================================================================================================
methylGSA_wrapper <- function(cpg.pval, array, method = "GSEA", GS, add_anno_for_sig = T, 
                              use_msigdbr = F, fun = "methylRRA", minsize = 100, maxsize = 500){
  
  if(add_anno_for_sig) {
    CpGs <- getAnnot(array.type = array) %>% as.data.frame()
  }
  
  results <- plyr::llply(
    GS,
    .fun = function(gs.type) {
      if(use_msigdbr) {
        pathway_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = toupper(gs.type))
        pathway <- unique(pathway_df$gs_name)
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
      
      fn <- eval(parse(text = fun))
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

  
      res_sig <- res %>% filter(.,padj < 0.05)
      
      # Add annotations
      if(add_anno_for_sig) {
        res_anno <- plyr::alply(
          res_sig,
          .margins = 1,
          .fun = function(pathway){

            if(method == "ORA") {
              pathway <- pathway %>% 
                separate_rows(overlap, sep = ",") 
              core <- "overlap"
            }
            if(method == "GSEA") {
              pathway <- pathway %>% 
                separate_rows(core_enrichment, sep = "/") 
              core <- "core_enrichment"
            }
           
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
        return(
          c(list(results = res),
            res_anno)
        )
       
      }
      
      res
     
    }
  )
  
  names(results) <- GS
  results

}

# Get annotation (source function from methylGSA)
getAnnot = function(array.type, group = "all"){
  if(array.type=="450K"){
    FullAnnot = tryCatch({
      getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    },
    error = function(e){
      stop("IlluminaHumanMethylation450kanno.ilmn12.hg19 needs to
be installed and loaded before running methylglm/methylRRA")
    })
  }else{
    FullAnnot = tryCatch({
      getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    },
    error = function(e){
      stop("IlluminaHumanMethylationEPICanno.ilm10b4.hg19 needs to
be installed and loaded before running methylglm/methylRRA")
    })
  }
  
  FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group")]
  FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
  FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]
  ## get the first gene in each USCS_RefGene_Name
  temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Name,split=";"),
                '[', 1, FUN.VALUE=character(1))
  FullAnnot$UCSC_RefGene_Name = temp
  ## get the first gene group in each UCSC_RefGene_Group
  temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Group,split=";"),
                '[', 1, FUN.VALUE=character(1))
  FullAnnot$UCSC_RefGene_Group = temp
  
  if(group == "body"){
    FullAnnot = 
      FullAnnot[FullAnnot$UCSC_RefGene_Group%in%c("Body", "1stExon"),]
  }
  
  if(group == "promoter1"){
    FullAnnot = FullAnnot[grepl("TSS",FullAnnot$UCSC_RefGene_Group),]
  }
  
  if(group == "promoter2"){
    FullAnnot = 
      FullAnnot[FullAnnot$UCSC_RefGene_Group%in%c("TSS200", "TSS1500", 
                                                  "1stExon", "5'UTR"),]
  }
  
  return(FullAnnot)
}