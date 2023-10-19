suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(stringr)
  library(corrplot)
  library(ComplexHeatmap)
  library(patchwork)
})

w.corr <- read_tsv(snakemake@input[['gene_loading_corr']])
rownames(w.corr) <- names(w.corr)

w.dist <- read_tsv(snakemake@input[['gene_loading_dist']])
rownames(w.dist) <- names(w.dist)

sig.num = length(w.corr) / 2
w.corr.mtx <- w.corr[1:sig.num, (sig.num+1):ncol(w.corr)] %>% as.matrix()
w.corr.mtx[w.corr.mtx < 0] = 0
sig.num = length(w.dist) / 2
w.dist.mtx <- w.dist[1:sig.num, (sig.num+1):ncol(w.dist)] %>% as.matrix()

validated.sig.df <- read_tsv(snakemake@input[['validated_sig_df']])

map.corr <- str_split(validated.sig.df[['validation.corr.1']], " ", simplify = T)[,2] %>% as.numeric()
map.dist <- str_split(validated.sig.df[['validation.dist.1']], " ", simplify = T)[,2] %>% as.numeric()

# plot gene loading correlations between discovery and validation
png(snakemake@output[['gene_loading_corr_plot']], width = 8, height = 6, units = "in", res = 300)
Heatmap((w.corr.mtx[,map.corr]), name = "Correlation",
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        row_labels = paste("discovery", seq_len(nrow(w.corr.mtx)), sep = " "),
        cluster_rows = F, cluster_columns = F)
dev.off()

# plot gene loading distance between discovery and validation
png(snakemake@output[['gene_loading_dist_plot']], width = 8, height = 6, units = "in", res = 300)
Heatmap((w.dist.mtx[,map.dist]), name = "Minkowski distance",
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        row_labels = paste("discovery", seq_len(nrow(w.dist.mtx)), sep = " "),
        cluster_rows = F, cluster_columns = F)
dev.off()


