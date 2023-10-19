suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(patchwork)
  library(tibble)
  library(scater)
  library(Nebulosa)
})

sce <- readRDS(snakemake@input[['sce']])
score_dframe <- readRDS(snakemake@input[['dframe']])

dimred = snakemake@params[['dimred_to_plot']]
sce_assay_to_plot = snakemake@params[['sce_assay_to_plot']]
singler_label_field = snakemake@params[['celltype_label_field']]
colour_label_field = snakemake@params[['metadata_to_color']]
individual_dimred_plot_dir = snakemake@params[['individual_celltype_score_dimred_dir']]

# plot annotated dimred plots, coloured by cell type labels and a metadata field
p1 <- plotReducedDim(sce,
                     dimred = dimred,
                     by_exprs_values = sce_assay_to_plot,
                     colour_by = singler_label_field) +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(colData(sce)[[singler_label_field]]))),
                     breaks = unique(colData(sce)[[singler_label_field]])) +
  labs(colour = singler_label_field)

p2 <- plotReducedDim(sce,
                     dimred = dimred,
                     by_exprs_values = sce_assay_to_plot,
                     colour_by = colour_label_field)

ggsave(snakemake@output[['dimred_cell_type']], (p1 + p2 + plot_layout(ncol = 2)), device = "png", width = 14, height = 7, units = "in", dpi = "retina")

# plot
celltypes <- score_dframe@listData$scores %>% colnames()

plist <- lapply(celltypes, function(celltype) {
  colData(sce)[['assignment score']] <- score_dframe@listData$scores[,celltype]
  plotReducedDim(sce,
                 dimred = dimred,
                 by_exprs_values = sce_assay_to_plot,
                 colour_by = 'assignment score') +
    labs(title = celltype)
})
names(plist) <- celltypes

if (!dir.exists(individual_dimred_plot_dir)) dir.create(individual_dimred_plot_dir, recursive = T)

lapply(celltypes, function(celltype) {
  ggsave(paste0(individual_dimred_plot_dir, dimred, "-celltype-assignment-score-", celltype, ".png"), (p1 + plist[[celltype]] + plot_layout(ncol = 2)), device = "png", width = 14, height = 7, units = "in", dpi = "retina")
})



