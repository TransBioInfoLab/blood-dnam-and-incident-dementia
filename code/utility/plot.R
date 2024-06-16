#######################################################################################################
# =================================================================================================== #
# Function for plot 
# =================================================================================================== #
#######################################################################################################
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(gtable)
library(grid)
library(cowplot)
library(meta)
# =================================================================================================== 
# Manhattan plot
# =================================================================================================== 
plot_manh <- function(results, annotated, colored, thres = 0.05, top_n = 10) {
  
  don <- results %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(pos)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(results, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each CpG
    arrange(chr, pos) %>%
    mutate( POScum=pos+tot) 

  annotated_df <- results %>% 
    filter(GREAT %in% annotated) %>% 
    dplyr::select(pVal.final, GREAT, CpG) %>% 
    group_by(GREAT) %>% 
    slice_min(pVal.final, with_ties = F)
  
  don <- don %>%
    # Add annotation information
    mutate( is_annotate=ifelse(CpG %in% annotated_df$CpG & fdr <= thres, "yes", "no"),
            is_red = ifelse(CpG %in% colored, "yes", "no")) 
  
  thres <- results %>%
    filter(fdr <= thres) %>%
    slice_max(pVal.final) %>%
    pull(pVal.final)
    
  # Prepare X axis
  axisdf <- don %>% 
    group_by(chr) %>% 
    summarize(center=( max(POScum) + min(POScum) ) / 2 ) %>%
    mutate(chr = ifelse(chr %in% c(11,13,15,17,19,21), "", chr))

  scaleFUN <- function(x) sprintf("%.0f", x)
  ggplot(don, aes(x=POScum, y=-log10(pVal.final))) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), size=1) +
    scale_color_manual(values = rep(c("grey20", "blue2"), 22 )) +
    geom_point(data = don %>% filter(is_red == "yes"), mapping = aes(x = POScum, y = -log10(pVal.final)), color = "red") +
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,-log10(1e-9)), labels=scaleFUN) +     # remove space between plot area and x axis
    
    # Add highlighted points
    geom_hline(aes(yintercept = -log10(thres)), color="red") +
    annotate("text", x = axisdf$center[2], y = -log10(1e-8), label = "FDR < 0.05") +
    annotate("segment", x = 0.6 * axisdf$center[1], xend = 1.4*axisdf$center[1], y = -log10(1e-8), yend = -log10(1e-8),
             colour = "red") + 
    # Add label using ggrepel to avoid overlapping
    geom_text_repel(data=subset(don, is_annotate=="yes"), 
                     aes(label=GREAT), 
                     size=4.5, fontface = "bold") +
    xlab("Chromosome") +
    ylab(expression("-log"["10"]~"(p)")) + 
    # Custom the theme:
    theme_classic() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    )
}
# =================================================================================================== 
# Barplot for pathway analysis
# =================================================================================================== 
plot_pathway <- function(results, thres_fdr = 0.05, title, ylim = 1e-5) {
  
  results <- results %>%
    filter(padj < thres_fdr) %>%
    mutate(logp = -log10(pvalue))
  
  
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  results <- results %>%
    mutate(ordering = pvalue,
           Description = fct_reorder(Description, ordering, .desc = F))
  cbf_1 <- c( "#E69F00", "#56B4E9", "#009E73", 
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  ggplot(results, aes(x = Description, y = logp)) + 
    geom_col(fill = "#0072B2", orientation = "x") +

    coord_flip() +
    scale_y_continuous(expand = c(0, 0), limits = c(0,-log10(ylim)), labels=scaleFUN) +

    theme_classic() +
    ylab(expression("-log"["10"]~"(p)")) + 
    xlab("") + 
    ggtitle(title) +
    theme(
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
      plot.title = element_text(size = 15, face = "bold"),
      strip.text.y = element_blank(),
      strip.placement = "outside",
      axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      legend.title = element_blank()
    )
  
  
  
}
# =================================================================================================== 
# Forest plot for meta analysis
# =================================================================================================== 
plot_forest <- function(results_list, cpg_order, effect, se, sm = 'HR', combine = T, ncol, nrow,  ...) {

  GREAT <- results_list[[1]]$GREAT_annotation

  # Match to the common cpg
  plot_list <- foreach::foreach(
    i = cpg_order,
    .errorhandling = "remove"
  ) %do% {

    dat <- purrr::map(results_list, ~.[.$cpg == i,])
    dat <- purrr::reduce(dat, rbind) 
    
    # Meta analysis
    f <- metagen(
      TE = dat[[effect]],
      seTE = dat[[se]],
      sm = sm,
      ...
    )

    grid.newpage()
    forest (
      f, 
      comb.fixed = TRUE, 
      comb.random = FALSE, 
      print.tau2 = TRUE,
      # sortvar = TE, # sorting variable
      digits.sd = 3,
      digits.se = 3,
      squaresize = 1, # size of the squares in the plot
      digits = 3,
      ddigits = 3, 
      digits.tau2 = 2,
      col.square = "blue",
      col.inside.fixed = "red",
      col.diamond = "red", 
      text.random = "Fixed effect model",
      leftcols = c("studlab"),
      rightcols = c("ci"),
      print.I2.ci = TRUE,
      smlab = i,
      colgap.forest.left = unit(4,"cm")
    )
    grid.grab()
  }

  if(combine) {
    ggpubr::ggarrange(
      plotlist = plot_list,
      ncol = ncol,
      nrow = nrow
    ) 
  } else {
    plot_list
  }
  
} 
# =================================================================================================== 
# Scatterplot for brain-blood correlation analysis
# =================================================================================================== 
plot_scatter <- function(x1, x2, data, ...) {
  ggpubr::ggscatter(
    data = data,
    x = x1,
    y = x2,
    cor.coef = T,
    add = "reg.line",
    add.params = list(color = "navy"),
    cor.method = "spearman",
    ...
  )
}
# =================================================================================================== 
# Dotplot for MRS analysis
# =================================================================================================== 
plot_dot_mrs <- function(data, CpGs = NA, td = NA){ 
  
  p <- ggpubr::ggdotchart(
    data,
    x = "Model",
    y = "AUC",
    facet.by = "Conversion",
    add.params = list(color = "lightgray", size = 1),
    font.label = list(color = "black", size = 8, vjust = 0.5, face = "bold"),     
    fill = "Age",
    sorting = "none",
    color = "Age",
    title = paste0("Training: ", td,  "\nCpGs: ", CpGs),
    label = round(data$AUC,digits = 3),
    add = "segments",                             # Add segments from y = 0 to dots
    rotate = TRUE,                               # Rotate vertically
    dot.size = 8,
    repel = T
  )  + ggsci::scale_fill_jco() + 
    ggsci::scale_color_jco() + 
    ylim(0.45,1) + 
    ggpubr::theme_cleveland()
  
  return(p)
}