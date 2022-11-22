suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
})

w.markers <- read_tsv(snakemake@input[['sig_known_gene_loading_mtx']])

# remove genes not in signature
w.markers.cleaned <- w.markers %>%
  drop_na(!split)

row_split <- w.markers$split
row_split_cleaned <- w.markers.cleaned$split

w.markers$split <- NULL
w.markers.cleaned$split <- NULL

gene.holder <- w.markers$gene
w.markers$gene <- NULL
rownames(w.markers) <- gene.holder

gene.holder <- w.markers.cleaned$gene
w.markers.cleaned$gene <- NULL
rownames(w.markers.cleaned) <- gene.holder

rm(gene.holder)

w.markers <- as.matrix(w.markers)
w.markers.cleaned <- as.matrix(w.markers.cleaned)

col_ha <- columnAnnotation(group = str_split(colnames(w.markers), pattern = " ", simplify = T)[,1])

heatmap_full <- Heatmap(w.markers, 
                        name = "Gene loading", 
                        cluster_rows = F, 
                        row_split = row_split, 
                        row_title_rot = 90,
                        row_title_gp = gpar(fontsize = 12),
                        bottom_annotation = col_ha)

heatmap_cleaned <- Heatmap(w.markers.cleaned, 
                           name = "Gene loading", 
                           cluster_rows = T, 
                           row_split = row_split_cleaned, 
                           row_title_rot = 90,
                           row_title_gp = gpar(fontsize = 12),
                           bottom_annotation = col_ha)

# save the plots
png(snakemake@output[['sig_known_gene_loading_plot_full']], width = 7, height = 18, units = "in", res = 300)
draw(heatmap_full,
     merge_legend = T)
dev.off()

png(snakemake@output[['sig_known_gene_loading_plot_cleaned']], width = 7, height = 18, units = "in", res = 300)
draw(heatmap_cleaned,
     merge_legend = T)
dev.off()

# =============================================================================

w.markers <- read_tsv(snakemake@input[['sig_known_gene_loading_mtx_validated']])

# remove genes not in signature
w.markers.cleaned <- w.markers %>%
  drop_na(!split)

row_split <- w.markers$split
row_split_cleaned <- w.markers.cleaned$split

w.markers$split <- NULL
w.markers.cleaned$split <- NULL

gene.holder <- w.markers$gene
w.markers$gene <- NULL
rownames(w.markers) <- gene.holder

gene.holder <- w.markers.cleaned$gene
w.markers.cleaned$gene <- NULL
rownames(w.markers.cleaned) <- gene.holder

rm(gene.holder)

w.markers <- as.matrix(w.markers)
w.markers.cleaned <- as.matrix(w.markers.cleaned)

col_ha <- columnAnnotation(group = str_split(colnames(w.markers), pattern = " ", simplify = T)[,1])

heatmap_full <- Heatmap(w.markers, 
                        name = "Gene loading", 
                        cluster_rows = F, 
                        row_split = row_split, 
                        row_title_rot = 90,
                        row_title_gp = gpar(fontsize = 12),
                        bottom_annotation = col_ha)

heatmap_cleaned <- Heatmap(w.markers.cleaned, 
                           name = "Gene loading", 
                           cluster_rows = T, 
                           row_split = row_split_cleaned, 
                           row_title_rot = 90,
                           row_title_gp = gpar(fontsize = 12),
                           bottom_annotation = col_ha)

# save the plots
png(snakemake@output[['sig_known_gene_loading_plot_full_validated']], width = 7, height = 18, units = "in", res = 300)
draw(heatmap_full,
     merge_legend = T)
dev.off()

png(snakemake@output[['sig_known_gene_loading_plot_cleaned_validated']], width = 7, height = 18, units = "in", res = 300)
draw(heatmap_cleaned,
     merge_legend = T)
dev.off()

