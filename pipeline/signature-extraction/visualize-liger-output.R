suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(rliger)
  library(sjstats)
  library(ggpubr)
  library(ggsci)
  library(tidyr)
  library(ComplexHeatmap)
  library(fastcluster)
})

liger <- readRDS(snakemake@input[['liger']])

# plot by cohort and cluster
all.plots <- plotByDatasetAndCluster(liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
ggsave(snakemake@output[['dataset_cluster_plot']], all.plots[[1]] + all.plots[[2]], width = 12, height = 7, units = "in")

liger <- runTSNE(liger)

# plot gene loadings 
if (!dir.exists(snakemake@params[['gene_loadings_plot_dir']])) dir.create(snakemake@params[['gene_loadings_plot_dir']], recursive = T)

datasets <- lapply(liger@raw.data, ncol) %>% unlist() %>% sort(decreasing = T) 
if (length(datasets) > 1) {
  dataset1 <- names(datasets)[1]
  dataset2 <- names(datasets)[2]
  
  gene_loadings <- plotGeneLoadings(liger, 
                                    dataset1 = dataset1,
                                    dataset2 = dataset2,
                                    do.spec.plot = FALSE, 
                                    return.plots = TRUE)
  
  lapply(seq(length(gene_loadings)), function(i) {
    ggsave(paste0(snakemake@params[['gene_loadings_plot_dir']], snakemake@wildcards[['subtype']], '-', snakemake@wildcards[['condition']], "-signature-", i, ".png"), gene_loadings[[i]], width = 7, height = 8, units = "in")
  })
} else {
  dataset1 <- names(datasets)[1]
  dataset2 <- names(datasets)[1]
}



# plot signature loadings
H.norm <- liger@H.norm
colnames(H.norm) <- paste("Signature", seq(ncol(H.norm)), sep = " ")

annotation_row <- data.frame(cohort = liger@cell.data$dataset)
rownames(annotation_row) <- rownames(H.norm)

annotation_col <- setNames(pal_jco()(length(unique(annotation_row$cohort))), unique(annotation_row$cohort))
col_fun = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

fh = function(x) fastcluster::hclust(dist(x))
p1 <- Heatmap(H.norm, 
              col = col_fun,
              show_row_names = F,
              name = "signature loading",
              right_annotation = rowAnnotation(df = annotation_row, 
                                               col = list(cohort = annotation_col)),
              row_split = annotation_row$cohort,
              cluster_rows = fh, 
              cluster_columns = fh)

# plot gene loading
W <- liger@W
rownames(W) <- paste("Signature", seq(nrow(W)), sep = " ")

p2 <- pheatmap(W, show_colnames = F, main = "gene loading")

# save the plots
png(snakemake@output[['signature_loading_plot']], width = 7, height = 12, units = "in", res = 300)
p1
dev.off()

png(snakemake@output[['gene_loading_plot']], width = 12, height = 7, units = "in", res = 300)
p2
dev.off()


