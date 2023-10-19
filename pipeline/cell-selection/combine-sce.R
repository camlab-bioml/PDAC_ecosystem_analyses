suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
})

source("utils/myhelpers.R")

sces <- lapply(snakemake@input[["sces"]], readRDS)

celltype_label_field <- snakemake@params[["celltype_label_field"]]
cell_type = snakemake@params[["cell_type"]]
print(paste0("combining: the ", snakemake@wildcards[["cohort"]], " cohort, for ", cell_type, " cells, using ", celltype_label_field, " labels."))
print("cell types used for combining:")
print(sapply(sces, function(sce) unique(colData(sce)[[celltype_label_field]])))

# Calculate the total number of cells per sample before combining
num_of_cells_before <- lapply(sces, function(sce) {
  if (ncol(sce) > 0) {
    df <- as.data.frame(table(sce$sample))
    df$condition <- unique(colData(sce)[[celltype_label_field]])
    df
  } else {
    data.frame(Var1 = character(), Freq = numeric(), condition = character())
  }
})
num_of_cells_before <- purrr::reduce(num_of_cells_before, rbind)

# combine to cell_type
sces <- lapply(sces, function(sce) {
  rowData(sce)[["seurat_variableFeatures_vst_varianceStandardized"]] <- NULL
  rowData(sce)[["seurat_variableFeatures_vst_mean"]] <- NULL
  reducedDim(sce, "Harmony") <- NULL
  sce
})
sce <- Reduce(cbind, sces)
colData(sce)[[celltype_label_field]] <- cell_type

# Calculate the total number of cells per sample after subsetting
if (ncol(sce) > 0) {
  num_of_cells_after <- as.data.frame(table(sce$sample))
  num_of_cells_after$condition <- paste(cell_type, "(combined)", sep = " ")
} else {
  num_of_cells_after <- data.frame(Var1 = character(), Freq = numeric(), condition = character())
}

# plot grouped bar plot to show number of cells in each sample before and after subsetting
num_of_cells <- bind_rows(num_of_cells_before, num_of_cells_after)
num_of_cells_plot <- ggplot(num_of_cells, aes(x = Var1, y = Freq, fill = condition)) +
  geom_bar(stat = 'identity', position = "dodge") +
  #facet_wrap(~group) +
  labs(x = 'Sample', y = 'Number of cells', title = paste0('Number of cells for all the samples before and after combining to ', cell_type, ' compartment')) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

png(snakemake@output[['cells_per_sample_plot']], width = 18, height = 8, units = 'in', res = 321)
num_of_cells_plot
dev.off()

#print(paste0(sum(num_of_cells_after$Freq) - sum(num_of_cells_before$Freq), " cells added by combining"))
print(paste0(sum(num_of_cells_after$Freq), " cells in total after combining"))

# magic that makes do_dimred (singleCellTK::runSeurat* functions) work
sce <- SingleCellExperiment(assays = assays(sce), rowData = rowData(sce), colData = colData(sce), reducedDims = reducedDims(sce), altExps = altExps(sce))

# if any cells were removed, redo PCA, Harmony, UMAP, and TSNE
if (sum(num_of_cells_before$Freq) > sum(num_of_cells_after$Freq) & sum(num_of_cells_after$Freq > 20)) {
  sce <- do_dimred(sce, harmonize = T, batch_col = "sample", n_cores = snakemake@threads)
}

saveRDS(sce, snakemake@output[['sce']])
names(num_of_cells) <- c("sample", "ncells", "condition")
write_tsv(num_of_cells, snakemake@output[['cells_per_sample']])


