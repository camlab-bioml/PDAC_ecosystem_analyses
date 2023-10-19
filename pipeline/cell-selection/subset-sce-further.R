suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
})

source("utils/myhelpers.R")

sce <- readRDS(snakemake@input[['sce']])

celltype_label_field <- snakemake@params[['celltype_label_field']]
subtype = snakemake@params[['subtype']]
print(paste0("subsetting: the ", snakemake@wildcards[['cohort']], " cohort, for ", subtype, " cells, using ", celltype_label_field, " labels."))
print("cell types available for subsetting:")
print(unique(colData(sce)[[celltype_label_field]]))

# Calculate the total number of cells per sample before subsetting
num_of_cells_before <- as.data.frame(table(sce$sample))
num_of_cells_before$condition <- paste0('Before subset to ', subtype)

# subset to immune cells
sce <- sce[,grepl(subtype, colData(sce)[[celltype_label_field]])]

# Calculate the total number of cells per sample after subsetting
if (ncol(sce) > 0) {
  num_of_cells_after <- as.data.frame(table(sce$sample))
  num_of_cells_after$condition <- paste0("After subset to ", cell_type)
} else {
  num_of_cells_after <- data.frame(Var1 = character(), Freq = numeric(), condition = character())
}

# plot grouped bar plot to show number of cells in each sample before and after subsetting
num_of_cells <- bind_rows(num_of_cells_before, num_of_cells_after)
num_of_cells_plot <- ggplot(num_of_cells, aes(x = Var1, y = Freq, fill = condition)) +
  geom_bar(stat = 'identity', position = "dodge") +
  #facet_wrap(~group) +
  labs(x = 'Sample', y = 'Number of cells', title = paste0('Number of cells for all the samples before and after subsetting to ', subtype, ' compartment')) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

png(snakemake@output[['cells_per_sample_plot']], width = 1800, height = 800)
num_of_cells_plot
dev.off()

print(paste0(sum(num_of_cells_before$Freq) - sum(num_of_cells_after$Freq), " cells removed by subsetting"))
print(paste0(sum(num_of_cells_after$Freq), " cells kept by subsetting"))

# magic that makes do_dimred (singleCellTK::runSeurat* functions) work
sce <- SingleCellExperiment(assays = assays(sce), rowData = rowData(sce), colData = colData(sce), reducedDims = reducedDims(sce), altExps = altExps(sce))

# if any cells were removed, redo PCA, Harmony, UMAP, and TSNE
if (sum(num_of_cells_before$Freq) > sum(num_of_cells_after$Freq) & sum(num_of_cells_after$Freq > 20)) {
  sce <- do_dimred(sce, harmonize = T, batch_col = "sample", n_cores = snakemake@threads)
}

saveRDS(sce, snakemake@output[['sce']])
names(num_of_cells) <- c("sample", "ncells", "condition")
write_tsv(num_of_cells, snakemake@output[['cells_per_sample']])


