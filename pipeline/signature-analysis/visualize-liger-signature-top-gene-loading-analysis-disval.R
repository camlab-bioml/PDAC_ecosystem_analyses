suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
})

w <- read_tsv(snakemake@input[['sig_top_gene_loading_mtx']])

ntop <- snakemake@params[["num_top_genes"]] %>% as.numeric()

gene.holder <- w$gene
w$gene <- NULL
rownames(w) <- gene.holder
rm(gene.holder)

group <- str_split(colnames(w), pattern = " ", simplify = T)[,1]
row_ha <- rowAnnotation(Group = group)
w.mtx <- as.matrix(w)

png(filename = snakemake@output[['sig_top_gene_loading_heatmap']], width = ntop*2, height = 7, units = "in", res = 321)
set.seed(42)
Heatmap(t(w.mtx), 
        name = "Gene loading",
        left_annotation = row_ha,
        show_column_names = FALSE,
        col = viridisLite::viridis(100, option = "C"))
dev.off()
