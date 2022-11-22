suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(magrittr)
  library(tidyverse)
  library(scater)
  library(Nebulosa)
})

seu <- readRDS(snakemake@input[['seu']])
sce <- readRDS(snakemake@input[['sce']])
score_mtx <- read_tsv(snakemake@input[['assignment_score']])

dimred = snakemake@params[['dimred_to_plot']]
annot_level = snakemake@params[['annot_level']]
sce_assay_to_plot = snakemake@params[['sce_assay_to_plot']]
mk = snakemake@params[['marker']]
individual_plots_dir = snakemake@params[['individual_plots_dir']]

output = snakemake@output[['expression_plot']]
print(output)

if(!dir.exists(individual_plots_dir)) dir.create(individual_plots_dir, recursive = T)

if (!(mk %in% rownames(sce))) {
  png(output, width = 1600,height = 1600)
  plot.new()
  dev.off()
} else {
  # normalizing expression to reference for better plotting
  seu <- NormalizeData(seu)
  Idents(seu) <- paste0("predicted.", annot_level)
  # plot marker expression on reduced dimensions
  p1 <- FeaturePlot(seu, reduction = dimred, features = mk)
  rm(seu)
  
  dimred = dimred %>% toupper()
  p2 <- plotReducedDim(sce, dimred, by_exprs_values = sce_assay_to_plot, colour_by = mk)
  
  if(max(logcounts(sce[mk,])) <= 0) {
    p3 <- p2
  } else {
    p3 <- plot_density(sce, mk, slot = sce_assay_to_plot, reduction = dimred, method = "ks", pal = "viridis")
  }
  
  # violin plots: marker expression against assigned cell types
  p4 <- plotExpression(sce, features = mk, 
                       x = paste0("predicted.", annot_level),
                       colour_by = paste0("predicted.", annot_level), 
                       exprs_values = sce_assay_to_plot) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = "Cell type assigned")
  
  # plot marker expression against predicted cell type scores
  cell_types <- score_mtx %>%
    select(!matches("cell_id|cell_type")) %>% 
    names()
  
  dflist <- lapply(cell_types, function(ct) {
    df <- data.frame(Expression = as.numeric(logcounts(sce[mk,])), 
                     Predicted_score = as.numeric(score_mtx[[ct]]),
                     Cell_type = ct)
  })
  df <- Reduce(bind_rows, dflist)
  
  # plot scatter plots with correlation between marker expression and predicted cell type scores
  p5 <- ggplot(df, aes(x = Expression, y = Predicted_score)) +
    geom_point(shape = 1) + 
    facet_wrap(vars(Cell_type)) + 
    ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), color = "red", geom = "label") +
    ggpubr::theme_pubr() + 
    labs(x = paste0(mk, " Expression"), y = "Predicted Score")
  
  # arrange plots for output
  png(output, width = 1600, height = 1600)
  print(p1 + p3 + p4 + p5 + plot_layout(ncol = 2))
  dev.off()
  
  # save individual plots
  png(paste0(individual_plots_dir, "expression.png"), width = 800, height = 800)
  print(p1)
  dev.off()
  
  png(paste0(individual_plots_dir, "expression-density.png"), width = 800, height = 800)
  print(p3)
  dev.off()
  
  png(paste0(individual_plots_dir, "expression-violin.png"), width = 800, height = 800)
  print(p4)
  dev.off()
  
  png(paste0(individual_plots_dir, "expression-vs-predicted-score.png"), width = 800, height = 800)
  print(p5)
  dev.off()
}




