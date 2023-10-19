suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(magrittr)
  library(tidyverse)
  library(scater)
  library(Nebulosa)
})

sce <- readRDS(snakemake@input[['sce']])
score_dframe <- readRDS(snakemake@input[['dframe']])

dimred = snakemake@params[['dimred_to_plot']]
celltype_label_field = snakemake@params[['celltype_label_field']]
sce_assay_to_plot = snakemake@params[['sce_assay_to_plot']]
mk = snakemake@params[['marker']]
individual_plots_dir = snakemake@params[['individual_plots_dir']]

output = snakemake@output[['expression_plot']]
print(output)

if(!dir.exists(individual_plots_dir)) dir.create(individual_plots_dir, recursive = T)

if (!(mk %in% rownames(sce))) {
  png(output, width = 18, height = 16, units = "in", res = 321)
  plot.new()
  dev.off()
} else {
  p1 <- plotReducedDim(sce, dimred, by_exprs_values = sce_assay_to_plot, colour_by = mk)
  
  if(max(logcounts(sce[mk,])) <= 0) {
    p2 <- p1
  } else {
    p2 <- plot_density(sce, mk, slot = sce_assay_to_plot, reduction = dimred, method = "ks", pal = "viridis")
  }
  
  # violin plots: marker expression against assigned cell types
  p3 <- plotExpression(sce, features = mk, 
                       x = celltype_label_field,
                       colour_by = celltype_label_field, 
                       exprs_values = sce_assay_to_plot) + 
    scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(colData(sce)[[celltype_label_field]]))),
                       breaks = unique(colData(sce)[[celltype_label_field]])) + 
    labs(colour = celltype_label_field) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = "Cell type assigned")
  
  # plot marker expression against predicted cell type scores
  score_mtx <- score_dframe@listData$scores %>% as.data.frame()
  cell_types <- names(score_mtx)
  
  dflist <- lapply(cell_types, function(ct) {
    df <- data.frame(Expression = as.numeric(logcounts(sce[mk,])), 
                     Predicted_score = as.numeric(score_mtx[[ct]]),
                     Cell_type = ct)
  })
  df <- Reduce(bind_rows, dflist)
  
  # plot scatter plots with correlation between marker expression and predicted cell type scores
  p4 <- ggplot(df, aes(x = Expression, y = Predicted_score)) +
    geom_point(shape = 1) + 
    facet_wrap(vars(Cell_type)) + 
    ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), color = "red", geom = "label") +
    ggpubr::theme_pubr() + 
    labs(x = paste0(mk, " Expression"), y = "Predicted Score")
  
  # arrange plots for output
  png(output, width = 18, height = 16, units = "in", res = 321)
  print(p1 + p2 + p3 + p4 + plot_layout(ncol = 2))
  dev.off()
  
  # save individual plots
  png(paste0(individual_plots_dir, "expression.png"), width = 9, height = 8, units = "in", res = 321)
  print(p1)
  dev.off()
  
  png(paste0(individual_plots_dir, "expression-density.png"), width = 9, height = 8, units = "in", res = 321)
  print(p2)
  dev.off()
  
  png(paste0(individual_plots_dir, "expression-violin.png"), width = 9, height = 8, units = "in", res = 321)
  print(p3)
  dev.off()
  
  png(paste0(individual_plots_dir, "expression-vs-predicted-score.png"), width = 9, height = 8, units = "in", res = 321)
  print(p4)
  dev.off()
}




