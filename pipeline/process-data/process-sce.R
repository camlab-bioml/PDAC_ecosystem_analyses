suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(tidyverse)
})

source("utils/myhelpers.R")

sce <- readRDS(snakemake@input[['sce']])
doublet_detection <- read_tsv(snakemake@input[['combined_doublet_detection']])

colData(sce)[["cohort"]] <- snakemake@wildcards[['cohort']]
print(snakemake@wildcards[['cohort']])

filtering_params <- snakemake@params[['filtering_params']]
print(filtering_params)

# Calculate the total number of cells per sample before filtering
num_of_cells_before <- as.data.frame(table(sce$sample))
num_of_cells_before$condition <- 'Before processing'

# find singlet ids
singlet_id <- doublet_detection %>%
  filter(classification == "Singlet") %>%
  pull(cell_id)

# remove doublets found by DoubletFinder
sce <- sce[,colnames(sce) %in% singlet_id]

# Calculate the total number of cells per sample after doublet removal
num_of_cells_singlet <- as.data.frame(table(sce$sample))
num_of_cells_singlet$condition <- 'After doublet removal'

# filter sce with custom QC filters
sce <- sce[, sce$subsets_mito_percent < filtering_params$mito_thresh & 
             sce$detected > filtering_params$detected_thresh & 
             sce$total > filtering_params$total_thresh]

# Calculate the total number of cells per sample after custom filtering
num_of_cells_after <- as.data.frame(table(sce$sample))
num_of_cells_after$condition <- 'After custom filters'

# plot grouped bar plot to show number of cells in each sample before and after processing
num_of_cells <- bind_rows(num_of_cells_before, num_of_cells_singlet, num_of_cells_after)
num_of_cells_plot <- ggplot(num_of_cells, aes(x = Var1, y = Freq, fill = condition)) +
  geom_bar(stat = 'identity', position = "dodge") +
  #facet_wrap(~group) +
  labs(x = 'Sample', y = 'Number of cells', title = 'Number of cells for all the samples before and after doublet removal and custom filtering') +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

png(snakemake@output[['cells_per_sample_plot']], width = 1800, height = 800)
num_of_cells_plot
dev.off()

print(paste0(sum(num_of_cells_before$Freq) - sum(num_of_cells_after$Freq), " cells removed by custom filters"))

# if any cells were removed, redo PCA, Harmony, UMAP, and TSNE
if (sum(num_of_cells_before$Freq) > sum(num_of_cells_after$Freq)) {
  sce <- do_dimred(sce, harmonize = T, batch_col = "sample", n_cores = snakemake@threads)
}

saveRDS(sce, snakemake@output[['sce']])
write_tsv(num_of_cells, snakemake@output[['cells_per_sample']])