suppressPackageStartupMessages({
  library(stringr)
  library(magrittr)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
})

source("utils/myhelpers.R")

scepath <- snakemake@input[['sce']]

cohorts <- str_split(string = scepath, pattern  = "-sce-", simplify = T)[,2]  
cohorts <- str_split(string = cohorts, pattern = ".rds", simplify = T)[,1]

print("cohorts available:")
print(cohorts)

cohorts_discovery = snakemake@params[['cohorts_discovery']]
cohorts_validation = snakemake@params[['cohorts_validation']]

print("discovery cohorts:")
print(cohorts_discovery)
print("validation cohorts:")
print(cohorts_validation)

# make sure the discovery-validation split makes sense for the full cohort list
stopifnot(all(c(cohorts_discovery, cohorts_validation) %in% cohorts))

# get the cohort sces
scelist <- lapply(scepath, readRDS)
names(scelist) <- cohorts

# store number of cells from each sample
num_of_cells <- lapply(scelist, function(sce) {
  holder <- as.data.frame(table(sce$sample))
  names(holder) = c("sample", "ncells")
  holder$cohort = sce$cohort[1]
  
  holder
})

num_of_cells <- Reduce(bind_rows, num_of_cells)
num_of_cells_discovery <- num_of_cells %>% filter(cohort %in% cohorts_discovery)
num_of_cells_validation <- num_of_cells %>% filter(cohort %in% cohorts_validation)

write_tsv(num_of_cells_discovery, snakemake@output[['cells_per_sample_discovery']])
write_tsv(num_of_cells_validation, snakemake@output[['cells_per_sample_validation']])

# group sces by discovery and validation, combine sces, and store discovery and validation group
scelist_discovery <- scelist[cohorts_discovery]
scelist_validation <- scelist[cohorts_validation]
scelist <- list(discovery = scelist_discovery,
                validation = scelist_validation)

# select dimred results to keep when cbind sces
dimred_to_keep = c("HARMONY", "PCA_ON_HARMONY", "UMAP_ON_HARMONY", "TSNE_ON_HARMONY", "INTEGRATED_DR", "REF.UMAP")

lapply(seq_along(scelist), function(i) {
  group = names(scelist)[[i]]
  sces = scelist[[i]]
  stopifnot(length(sces) >= 1)
  
  sces <- lapply(sces, function(sce) {
    rownames(sce) <- paste(rownames(sce), rowData(sce)[['ensembl_id']], sep = "_")
    sce
  })
  saveRDS(sces, snakemake@output[[paste0('scelist_', group)]])
  
  if (length(sces) > 1) {
    geneslist <- lapply(sces, rownames)
    common_genes <- Reduce(intersect, geneslist)
    sces <- lapply(sces, function(sce) {sce[common_genes,]})
  }
  
  sces <- lapply(sces, function(sce) {
    reducedDims(sce) <- reducedDims(sce)[reducedDimNames(sce) %in% dimred_to_keep]
    sce
  })
  sce <- Reduce(cbind, sces)
  rownames(sce) <- rowData(sce)[['Symbol']]
  
  # redo PCA, Harmony, UMAP, and TSNE; correcting for cohort level batch effect
  if (length(sces) > 1 & ncol(sce) > 20) {sce <- do_dimred(sce, harmonize = T, batch_col = "cohort", n_cores = snakemake@threads)}
  
  saveRDS(sce, snakemake@output[[paste0('sce_', group)]])

})


