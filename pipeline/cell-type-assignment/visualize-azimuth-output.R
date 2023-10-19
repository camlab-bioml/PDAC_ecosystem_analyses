suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(magrittr)
  library(tidyverse)
  library(ComplexHeatmap)
  library(fastcluster)
})

seu <- readRDS(snakemake@input[['seu']])
score_mtx <- read_tsv(snakemake@input[['assignment_score']])

dimred = snakemake@params[['dimred_to_plot']]

print("Dimensionality reductions available for plotting:")
Reductions(seu)
print("Dimensionality reduction used for plotting:")
print(dimred)

# plot projection of the data on the reference, coloured by samples and predicted cell types
p1 <- DimPlot(seu, reduction = dimred, group.by = paste0("predicted.", snakemake@params[['annot_level']]), label = TRUE, repel = T, label.size = 6) + NoLegend()
p2 <- DimPlot(seu, reduction = dimred, group.by = snakemake@params[['metadata_to_color']])

png(snakemake@output[['umap_cell_type']], width = 18, height = 8, units = "in", res = 321)
p1 + p2
dev.off()

# normalize seurat data for better plotting 
seu <- NormalizeData(seu)
Idents(seu) <- paste0("predicted.", snakemake@params[['annot_level']])

# peek at cell types predicted by Azimuth
cell_types <- seu@assays[[paste0("prediction.score.", snakemake@params[['annot_level']])]]@data %>% rownames()
print("cell types predicted by Azimuth: ")
print(cell_types)

# plot predicted cell type scores on reduced dimension
plist <- lapply(cell_types, function(ct) {
  p <- FeaturePlot(seu, 
              reduction = dimred,
              features = paste0("predictionscore", gsub(".", "", snakemake@params[['annot_level']], fixed = T), "_", ct)) +
    scale_color_gradient(limits = c(0,1), low = "grey80", high = "blue", name = "predicted score") + 
    labs(title = ct)
  p
})
names(plist) <- cell_types

lapply(names(plist), function(ct) {
  if (!dir.exists(snakemake@params[['individual_celltype_prediction_score_UMAP_dir']])) {
    dir.create(snakemake@params[['individual_celltype_prediction_score_UMAP_dir']], recursive = T)
  }
  png(paste0(snakemake@params[['individual_celltype_prediction_score_UMAP_dir']], 'UMAP-celltype-prediction-score-', ct,'.png'), 
      width = 12, height = 8, units = "in", res = 321)
  print(plist[[ct]])
  dev.off()
})

# save combined cell type prediction sore plot
p3 <- wrap_plots(plist, guides = "collect")
png(snakemake@output[['umap_prediction_score']], width = 18, height = 15, units = "in", res = 321)
p3
dev.off()

# plot prediction scores for single cells as heatmap
if(nrow(score_mtx) > 50000) {
  score_mtx <- sample_n(score_mtx, 40000)
}

mtx <- score_mtx %>% 
  select(!matches("cell_id|cell_type")) %>% 
  as.matrix() %>%
  t()
colnames(mtx) <- score_mtx[['cell_id']]

top_ha <- HeatmapAnnotation(Assigned = score_mtx[['cell_type']], 
                            name = "Assigned cell type", 
                            col = list(Assigned = structure(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(score_mtx[['cell_type']]))), 
                                                            names = unique(score_mtx[['cell_type']]))))
hm <- Heatmap(mtx, name = "Predicted score", 
              col = colorRampPalette(viridisLite::viridis(100))(100),
              show_column_names = F,
              top_annotation = top_ha,
              heatmap_legend_param = list(direction = "vertical"),
              column_split = score_mtx[['cell_type']],
              use_raster = T)

png(snakemake@output[['heatmap_prediction_score']], width = 12, height = 5, units = "in", res = 321)
draw(hm, merge_legend = T)
dev.off()



