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
w.dist <- as.matrix(w.dist)

# plot gene loading correlations between discovery and validation
png(snakemake@output[['gene_loading_corr_plot']], width = 7, height = 7, units = "in", res = 300)
corrplot(corr = as.matrix(w.corr))
dev.off()

# plot gene loading distance between discovery and validation
w.dist[w.dist == 0] <- NA
group <- str_split(colnames(w.dist), pattern = " ", simplify = T)[,1]

png(snakemake@output[['gene_loading_dist_plot']], width = 8, height = 6, units = "in", res = 300)
Heatmap(w.dist, name = paste(snakemake@params[['dist_method']], "distance", sep = " "),
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        cluster_rows = F, cluster_columns = F,
        row_split = group,
        column_split = group)
dev.off()

#sig.num = ncol(w.dist) / 2
#w.dist.mtx <- w.dist[(sig.num+1):nrow(w.dist),1:sig.num]
#pheatmap((w.dist.mtx), cluster_rows = F, cluster_cols = F)
