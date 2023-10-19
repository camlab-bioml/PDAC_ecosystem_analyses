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

w.corr <- read_tsv(snakemake@input[["validated_gene_loading_corr"]])
rownames(w.corr) <- names(w.corr)

z.corr <- read_tsv(snakemake@input[["validated_sig_loading_corr"]])
rownames(z.corr) <- names(z.corr)

w.corr.mtx <- w.corr %>% as.matrix()
z.corr.mtx <- z.corr %>% as.matrix()

# plot gene loading correlations between discovery and validation
png(snakemake@output[["validated_gene_loading_corr_plot"]], width = 8, height = 6, units = "in", res = 300)
Heatmap(w.corr.mtx, name = "Gene loading corr.",
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        cluster_rows = TRUE, cluster_columns = TRUE)
dev.off()

# plot gene loading distance between discovery and validation
png(snakemake@output[["validated_sig_loading_corr_plot"]], width = 8, height = 6, units = "in", res = 300)
Heatmap(z.corr.mtx, name = "Signature loading corr.",
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        cluster_rows = TRUE, cluster_columns = TRUE)
dev.off()


