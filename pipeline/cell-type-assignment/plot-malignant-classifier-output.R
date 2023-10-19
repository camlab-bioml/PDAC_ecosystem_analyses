suppressPackageStartupMessages({
        library(magrittr)
        library(tidyverse)
        library(sjstats)
        library(ggpubr)
        library(ggsci)
        library(stringr)
        library(ComplexHeatmap)
        library(Matrix)
        library(here)
        library(scales)
        library(gridExtra)
        library(bluster)
})

set.seed(123L)

# load tumor original expression sce
sce.tumor.results.list <- readRDS(snakemake@input[["sce_tumor_results_list"]])
sce.tumor.infercnv.list <- readRDS(snakemake@input[["sce_tumor_infercnv_list"]])

the.genes <- c("KRT19", "CFTR")

# cluster cells using infercnv smoothed expression
lapply(sce.tumor.infercnv.list, function(sce.tumor.infercnv) {
  p0 <- ggplot(data = as.data.frame(colData(sce.tumor.infercnv)) |> select(label, infer_cnv_var), aes(x = label, y = infer_cnv_var, color = label)) + 
    geom_boxplot() +
    scale_y_log10() +
    #scale_color_npg() +
    geom_hline(yintercept = mean(colData(sce.tumor.infercnv)$infer_cnv_var), linetype = 2) +
    #geom_pwc(aes(group = label), tip.length = 0, method = "t_test", label = "p.adj.format") +
    #stat_compare_means() + 
    theme_pubr()
  
  p <- plotReducedDim(sce.tumor.infercnv, "UMAP", colour_by="infer_cnv_var")
  p <- p + scale_colour_continuous(name = "infer_cnv_var", type = "viridis", 
                                   limits=c(0, 0.002), 
                                   oob=squish)
  gridExtra::grid.arrange(
    p0,
    p,
    plotReducedDim(sce.tumor.infercnv, "UMAP", colour_by="label"),
    ncol=3
  )
})

# plot KRT19 and CFTR expression 
the.genes = c("KRT19", "CFTR")

lapply(sce.tumor.infercnv.list, function (sce.tumor) {
  sce.tumor$greater_than_mean_infer_cnv_var <- sce.tumor$infer_cnv_var > mean(sce.tumor$infer_cnv_var)
  plotExpression(
    sce.tumor,
    features = the.genes,
    x = "greater_than_mean_infer_cnv_var",
    exprs_values = "logcounts",
    colour_by = "label",
    xlab = "Clusters",
    one_facet = TRUE,
    ncol = 1,
    scales = "free"
  )
})

# plots to help setting decision threshold
pdf(snakemake@output[['cell_cluster_plots']], width = 15, height = 15, bg = "white")
lapply(sce.tumor.infercnv.list, function(sce.tumor.infercnv) {
  sce.tumor <- sce.tumor.results.list[[(sce.tumor.infercnv$sample)[1]]]
  sce.tumor$label <- sce.tumor.infercnv$label
  sce.tumor$cnv_state_var <- sce.tumor.infercnv$cnv_state_var
  sce.tumor$greater_than_mean_infer_cnv_var <- sce.tumor$infer_cnv_var > mean(sce.tumor$infer_cnv_var, na.rm = T)
  sce.tumor$greater_than_mean_cnv_state_var <- sce.tumor$cnv_state_var > mean(sce.tumor$cnv_state_var, na.rm = T)
  sce.tumor$greater_than_1_cnv_state_var <- sce.tumor$cnv_state_var > 1
  
  p33 <- plotExpression(
    sce.tumor,
    features = the.genes,
    x = "label",
    exprs_values = "logcounts",
    colour_by = "cnv_state_var",
    xlab = "Clusters",
    one_facet = TRUE,
    ncol = 1,
    scales = "free"
  )
  
  p00 <- plotExpression(
    sce.tumor,
    features = the.genes,
    x = "greater_than_mean_infer_cnv_var",
    exprs_values = "logcounts",
    colour_by = "label",
    xlab = "greater_than_mean_infer_cnv_var",
    one_facet = TRUE,
    ncol = 1,
    scales = "free"
  )
  
  p11 <- plotExpression(
    sce.tumor,
    features = the.genes,
    x = "greater_than_mean_cnv_state_var",
    exprs_values = "logcounts",
    colour_by = "label",
    xlab = "greater_than_mean_cnv_state_var",
    one_facet = TRUE,
    ncol = 1,
    scales = "free"
  )
  
  p22 <- plotExpression(
    sce.tumor,
    features = the.genes,
    x = "greater_than_1_cnv_state_var",
    exprs_values = "logcounts",
    colour_by = "label",
    xlab = "greater_than_1_cnv_state_var",
    one_facet = TRUE,
    ncol = 1,
    scales = "free"
  )
  
  p1 <- ggplot(data = as.data.frame(colData(sce.tumor.infercnv)) |> select(label, cnv_state_var), aes(x = label, y = cnv_state_var, color = label)) + 
    geom_boxplot() +
    #scale_color_npg() +
    geom_hline(yintercept = mean(colData(sce.tumor.infercnv)$cnv_state_var, na.rm = T), linetype = 2) +
    geom_pwc(aes(group = label), tip.length = 0, method = "t_test", label = "p.adj.format") +
    #stat_compare_means() + 
    labs(title = (sce.tumor.infercnv$sample)[1]) +
    theme_pubr()
  
  p0 <- ggplot(data = as.data.frame(colData(sce.tumor.infercnv)) |> select(label, infer_cnv_var), aes(x = label, y = infer_cnv_var, color = label)) + 
    geom_boxplot() +
    #scale_color_npg() +
    geom_hline(yintercept = mean(colData(sce.tumor.infercnv)$infer_cnv_var, na.rm = T), linetype = 2) +
    geom_pwc(aes(group = label), tip.length = 0, method = "t_test", label = "p.adj.format") +
    #stat_compare_means() + 
    labs(title = (sce.tumor.infercnv$sample)[1]) +
    theme_pubr()
  
  p111 <- plotReducedDim(sce.tumor.infercnv, "UMAP", colour_by="cnv_state_var")
  p111 <- p111 + scale_colour_continuous(name = "cnv_state_var", type = "viridis", 
                                         #limits=c(0, 0.002), 
                                         oob=squish)
  
  p000 <- plotReducedDim(sce.tumor.infercnv, "UMAP", colour_by="infer_cnv_var")
  p000 <- p000 + scale_colour_continuous(name = "infer_cnv_var", type = "viridis", 
                                         #limits=c(0, 0.002), 
                                         oob=squish)
 tryCatch({
   gridExtra::grid.arrange(
    p000,
    p00,
    p0,
    p111,
    p11,
    p1,
    plotReducedDim(sce.tumor.infercnv, "UMAP", colour_by="label"),
    p22,
    p33,
    ncol=3
  )
 }, error = function(e) {
   gridExtra::grid.arrange(
    p000,
    p00,
    p0,
    p111,
    p11,
    p1,
    plotReducedDim(sce.tumor.infercnv, "UMAP", colour_by="label"),
    p33,
    ncol=3
  )
 })
})
dev.off()

# plot the expression of the genes in the classified cells
pdf(snakemake@output[['cell_marker_exprs_plots']], width = 10, height = 7, bg = "white")
lapply(sce.tumor.results.list, function(sce.tumor.results) {
  sce.tumor.results$cell_type <- sce.tumor.infercnv.list[[(sce.tumor.results$sample)[1]]]$cell_type
  sce.tumor.results$cell_type_cell_by_cell <- sce.tumor.infercnv.list[[(sce.tumor.results$sample)[1]]]$cell_type_cell_by_cell
  sce.tumor.results$label <- sce.tumor.infercnv.list[[(sce.tumor.results$sample)[1]]]$label
  p1 <- plotExpression(
    sce.tumor.results,
    features = the.genes,
    x = "cell_type",
    exprs_values = "logcounts",
    colour_by = "label",
    xlab = "cell_type",
    one_facet = TRUE,
    ncol = 1,
    scales = "free"
  ) + labs(title = (sce.tumor.results$sample)[1])
  p2 <- plotExpression(
    sce.tumor.results,
    features = the.genes,
    x = "cell_type_cell_by_cell",
    exprs_values = "logcounts",
    colour_by = "label",
    xlab = "cell_type by thresholding per cell",
    one_facet = TRUE,
    ncol = 1,
    scales = "free"
  ) + labs(title = (sce.tumor.results$sample)[1])

  gridExtra::grid.arrange(
    p1,
    p2,
    ncol = 2
  )
})
dev.off()