suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
})

w <- read_tsv(snakemake@input[['sig_top_gene_loading_mtx']])

gene.holder <- w$gene
w$gene <- NULL
rownames(w) <- gene.holder
rm(gene.holder)

group <- str_split(colnames(w), pattern = " ", simplify = T)[,1]
row_ha <- rowAnnotation(Group = group)
w.mtx <- as.matrix(w)

png(filename = snakemake@output[['sig_top_gene_loading_heatmap']], width = 15, height = 7, units = "in", res = 300)
Heatmap(t(w.mtx), 
        name = "Gene loading",
        left_annotation = row_ha)
dev.off()

# =============================================================================

w <- read_tsv(snakemake@input[['sig_top_gene_loading_mtx_validated']])

gene.holder <- w$gene
w$gene <- NULL
rownames(w) <- gene.holder
rm(gene.holder)

group <- str_split(colnames(w), pattern = " ", simplify = T)[,1]
row_ha <- rowAnnotation(Group = group)
w.mtx <- as.matrix(w)

png(filename = snakemake@output[['sig_top_gene_loading_heatmap_validated']], width = 15, height = 7, units = "in", res = 300)
Heatmap(t(w.mtx), 
        name = "Gene loading",
        left_annotation = row_ha)
dev.off()
