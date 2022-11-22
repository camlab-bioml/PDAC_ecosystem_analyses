suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(magrittr)
})

print(snakemake@params[['doublet_detection_dir']])

sce <- readRDS(snakemake@input[['sce']])
samples <- unique(sce$sample)
rm(sce)

doublet_detection_filenames <- paste0("DoubletFinder-", samples, ".tsv")

detections_per_sample <- lapply(doublet_detection_filenames, function(fname) {
  df <- read_tsv(paste0(snakemake@params[['doublet_detection_dir']], fname)) %>% 
    select(cell_id, contains("classifications")) %>%
    select(cell_id, last_col())
  names(df) <- c("cell_id", "classification")
  df
})

detections <- Reduce(bind_rows, detections_per_sample)

write_tsv(detections, snakemake@output[['combined_doublet_detection']])