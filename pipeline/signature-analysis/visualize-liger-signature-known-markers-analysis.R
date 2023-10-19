suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
  library(GeneOverlap)
})

w.markers <- read_tsv(snakemake@input[['sig_known_gene_loading_mtx']])

# remove genes not in signature
w.markers.cleaned <- w.markers %>%
  drop_na(!split)
w.markers.cleaned <- w.markers.cleaned[rowSums(w.markers.cleaned %>% select(where(is.numeric))) > 0, ]

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
w.markers.cleaned <- as.matrix(w.markers.cleaned) %>% t() %>% scale() %>% t()

col_ha <- columnAnnotation(group = str_split(colnames(w.markers), pattern = " ", simplify = T)[,1])

heatmap_full <- Heatmap(w.markers, 
                        name = "Gene loading", 
                        cluster_rows = FALSE, 
                        row_split = row_split, 
                        row_title_rot = 0,
                        row_title_gp = gpar(fontsize = 14, fontface = "bold"),
                        show_row_names = FALSE,
                        column_labels = gsub(" Rep| RepVal", "", colnames(w.markers)),
                        column_names_gp = gpar(fontsize = 15, fontface = "bold"),
                        col = viridisLite::viridis(100, option = "C"),
                        height = unit(13, "in")
                        #bottom_annotation = col_ha
                        )

heatmap_cleaned <- Heatmap(w.markers.cleaned, 
                           name = "Gene loading", 
                           cluster_rows = TRUE, 
                           row_split = row_split_cleaned, 
                           row_title_rot = 0,
                           row_title_gp = gpar(fontsize = 14, fontface = "bold"),
                           show_row_names = FALSE,
                           column_labels = gsub(" Rep| RepVal", "", colnames(w.markers.cleaned)),
                           column_names_gp = gpar(fontsize = 15, fontface = "bold"),
                           col = viridisLite::viridis(100, option = "C"),
                           height = unit(13, "in")
                           #bottom_annotation = col_ha
                           )

# save the plots
png(snakemake@output[['sig_known_gene_loading_plot_full']], width = 8, height = 18, units = "in", res = 360)
draw(heatmap_full,
     merge_legend = T)
dev.off()

png(snakemake@output[['sig_known_gene_loading_plot_cleaned']], width = 8, height = 18, units = "in", res = 360)
draw(heatmap_cleaned,
     merge_legend = T)
dev.off()


# plot GOM (heatmaps) ----------------------------------------------------------
gom.obj <- readRDS(snakemake@input[['sig_gom']])
pval.cutoff <- as.numeric(snakemake@params[['pvalue_cutoff_for_overlap_test']])

png(snakemake@output[['sig_gom_plot_oddsratio']], width = 22, height = 12, units = "in", res = 360)
tryCatch({
  drawHeatmap(gom.obj, 
              what = "odds.ratio", 
              log.scale = T, 
              adj.p = T, 
              cutoff = pval.cutoff, 
              ncolused = 5, 
              grid.col = "Blues", 
              note.col = "black")
}, error = function(e) {
  par(mar=c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, 
       paste("drawHeatmap() encountered this error:\n", e), 
       cex = 1, col = "black", family = "serif", font = 2, adj = 0.5)
})
dev.off()

png(snakemake@output[['sig_gom_plot_jaccardindex']], width = 22, height = 12, units = "in", res = 360)
tryCatch({
  drawHeatmap(gom.obj, 
              what = "Jaccard", 
              log.scale = T, # does not work for Jaccard index
              adj.p = T, 
              cutoff = pval.cutoff, 
              ncolused = 5, 
              grid.col = "Oranges", 
              note.col = "black")
}, error = function(e) {
  par(mar=c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, 
       paste("drawHeatmap() encountered this error:\n", e), 
       cex = 1, col = "black", family = "serif", font = 2, adj = 0.5)
})
dev.off()



