suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
})

source("utils/myhelpers.R")

sce <- readRDS(snakemake@input[['sce']])

annot_level <- snakemake@params[['annot_level']]
cell_type = snakemake@params[['cell_type']]
print(paste0("subsetting: the ", snakemake@wildcards[['cohort']], " cohort, for ", cell_type, " cells, using ", annot_level, " labels."))
print("cell types available for subsetting:")
print(unique(colData(sce)[[paste0("predicted.", annot_level)]]))

# Calculate the total number of cells per sample before subsetting
num_of_cells_before <- as.data.frame(table(sce$sample))
num_of_cells_before$condition <- paste0('Before subset to ', cell_type)

# subset to immune cells
sce <- sce[,grepl(cell_type, colData(sce)[[paste0("predicted.", annot_level)]])]

# Calculate the total number of cells per sample after subsetting
num_of_cells_after <- as.data.frame(table(sce$sample))
num_of_cells_after$condition <- paste0('After subset to ', cell_type)

# plot grouped bar plot to show number of cells in each sample before and after subsetting
num_of_cells <- bind_rows(num_of_cells_before, num_of_cells_after)
num_of_cells_plot <- ggplot(num_of_cells, aes(x = Var1, y = Freq, fill = condition)) +
  geom_bar(stat = 'identity', position = "dodge") +
  #facet_wrap(~group) +
  labs(x = 'Sample', y = 'Number of cells', title = paste0('Number of cells for all the samples before and after subsetting to ', cell_type, ' compartment')) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

png(snakemake@output[['cells_per_sample_plot']], width = 1800, height = 800)
num_of_cells_plot
dev.off()

print(paste0(sum(num_of_cells_before$Freq) - sum(num_of_cells_after$Freq), " cells removed by subsetting"))
print(paste0(sum(num_of_cells_after$Freq), " cells kept by subsetting"))

# if any cells were removed, redo PCA, Harmony, UMAP, and TSNE
if (sum(num_of_cells_before$Freq) > sum(num_of_cells_after$Freq) & sum(num_of_cells_after$Freq > 20)) {
  sce <- do_dimred(sce, harmonize = T, batch_col = "sample", n_cores = snakemake@threads)
}

saveRDS(sce, snakemake@output[['sce']])
names(num_of_cells) <- c("sample", "ncells", "condition")
write_tsv(num_of_cells, snakemake@output[['cells_per_sample']])


