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
  library(circlize)
  library(scales)
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

w.corr.mtx <- w.corr.mtx[,map.corr]
w.dist.mtx <- w.dist.mtx[,map.corr]

w.corr.mtx <- matrix(as.numeric(w.corr.mtx), ncol = ncol(w.corr.mtx), 
                     dimnames = list(paste("discovery", seq_len(nrow(w.corr.mtx)), sep = " "), 
                                     colnames(w.corr.mtx)))
w.dist.mtx <- matrix(as.numeric(w.dist.mtx) %>% rescale(), ncol = ncol(w.dist.mtx), 
                     dimnames = list(paste("discovery", seq_len(nrow(w.dist.mtx)), sep = " "),
                                     colnames(w.dist.mtx)))

# plot gene loading correlations between discovery and validation
ht.coor <- Heatmap(w.corr.mtx, name = "Correlation", rect_gp = gpar(type = "none"),
#col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        col = colorRamp2(seq(from = min(w.corr.mtx), to = max(w.corr.mtx), length.out = 100), viridisLite::viridis(100, option = "D")),
        cluster_rows = F, cluster_columns = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (i >= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "white"))
          }
        })

png(snakemake@output[['gene_loading_corr_plot']], width = 6, height = 4, units = "in", res = 360)
draw(ht.coor)
dev.off()

# plot gene loading distance between discovery and validation
ht.dist <- Heatmap(w.dist.mtx, name = "Minkowski distance", rect_gp = gpar(type = "none"),
        col = colorRamp2(seq(from = min(w.dist.mtx), to = max(w.dist.mtx), length.out = 100), rev(viridisLite::viridis(100, option = "E"))),
        cluster_rows = F, cluster_columns = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (i <= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "white"))
          }
        },
        show_column_names = FALSE)

png(snakemake@output[['gene_loading_dist_plot']], width = 6, height = 4, units = "in", res = 360)
draw(ht.dist)
dev.off()

# plot both gene loading correlations and distance between discovery and validation
png(snakemake@output[['gene_loading_corr_and_dist_plot']], width = 8, height = 6, units = "in", res = 360)
draw(ht.coor + ht.dist, ht_gap = unit(-80, "mm"))
dev.off()
