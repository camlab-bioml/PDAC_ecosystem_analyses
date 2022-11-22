suppressPackageStartupMessages({
  library(magrittr)
  library(SingleCellExperiment)
  library(scuttle)
  library(scater)
  library(SingleR)
  library(ggrepel)
  library(ggpubr)
})

pred <- readRDS(snakemake@input[['dframe']])
sce <- readRDS(snakemake@input[['sce']])

ncol_violin = snakemake@params[['ncol_violin']]
celltype_label_field = snakemake@params[['celltype_label_field']]

# plot SingleR predicted score for each cell type label
png(snakemake@output[['heatmap_assignment_score']], width = 800, height = 1200)
plotScoreHeatmap(pred, 
                 show.pruned = T, 
                 show.labels = T)
dev.off()

# plot distribution of delta between assigned label score and median score for single cells
png(snakemake@output[['violin_delta_distribution']], width = 1200, height = 800)
plotDeltaDistribution(pred,
                      show = "delta.med",
                      ncol = ncol_violin)
dev.off()

# plot distribution of scores for each label, grouped by label status
png(snakemake@output[['violin_score_distribution']], width = 1200, height = 800)
plotScoreDistribution(pred,
                      ncol = ncol_violin)
dev.off()

# plot marker expression heatmap for each cell type
all.markers <- metadata(pred)$de.genes

colormap <- structure(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(colData(sce)[[celltype_label_field]]))), 
                    names = unique(colData(sce)[[celltype_label_field]]))
colormap <- list(colormap)
names(colormap) <- celltype_label_field

## cell-related markers
plist <- lapply(names(all.markers), function(cell_type) {
  plotHeatmap(sce,
              color_columns_by = celltype_label_field,
              column_annotation_colors = colormap,
              order_columns_by = celltype_label_field,
              features = unique(unlist(all.markers[[cell_type]])),
              show_colnames = F)
})
names(plist) <- names(all.markers)

lapply(names(plist), function(ct) {
  if (!dir.exists(snakemake@params[['individual_celltype_marker_expression_heatmap_dir']])) {
    dir.create(snakemake@params[['individual_celltype_marker_expression_heatmap_dir']], recursive = T)
  }
  png(paste0(snakemake@params[['individual_celltype_marker_expression_heatmap_dir']], 'heatmap-celltype-marker-expression-', ct,'.png'), 
      width = 800, height = 2400)
  print(plist[[ct]])
  dev.off()
})


