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
num_of_cells <- list(discovery = num_of_cells_discovery,
                     validation = num_of_cells_validation)

# group sces by discovery and validation, combine sces, and store discovery and validation group
scelist_discovery <- scelist[cohorts_discovery]
scelist_validation <- scelist[cohorts_validation]
scelist <- list(discovery = scelist_discovery,
                validation = scelist_validation)

# select dimred results to keep when cbind sces
dimred_to_keep = c("UMAP_ON_HARMONY", "TSNE_ON_HARMONY", "INTEGRATED_DR", "REF.UMAP")

# the grouping of cohorts
lapply(seq_along(scelist), function(i) {
  group = names(scelist)[[i]]
  sces = scelist[[group]]
  sample_cell_counts = num_of_cells[[group]]
  
  # filter out samples with too few cells
  sample_cell_counts <- sample_cell_counts %>% 
    filter(ncells >= as.numeric(snakemake@params[['sample_cell_count_thres']]))
  
  # handle the sces
  sces <- lapply(sces, function(sce) {
    samples_to_keep <- sample_cell_counts %>% filter(cohort == sce[["cohort"]][1])
    samples_to_keep <- samples_to_keep[["sample"]]
    
    sce[,sce$sample %in% samples_to_keep]
  })
  
  for (c in names(sces)) {
    if (ncol(sces[[c]]) == 0) sces[[c]] <- NULL
  }
  rm(c)
  
  # make sure there is at least one cohort in the group
  stopifnot(length(sces) >= 1)
  
  # update gene names
  sces <- lapply(sces, function(sce) {
    rownames(sce) <- paste(rownames(sce), rowData(sce)[['ensembl_id']], sep = "_")
    sce
  })
  
  # get common genes and cell metadata DFrame across multiple cohorts
  if (length(sces) > 1) {
    geneslist <- lapply(sces, rownames)
    common_genes <- Reduce(intersect, geneslist)
    sces <- lapply(sces, function(sce) {sce[common_genes,]})
    
    coldata_fields_list <- lapply(sces, function(sce) names(colData(sce)))
    common_coldata_fields <- Reduce(intersect, coldata_fields_list)
    sces <- lapply(sces, function(sce) {
      for (field in setdiff(names(colData(sce)), common_coldata_fields)) {
        colData(sce)[[field]] <- NULL
      }
      sce
    })
  }
  
  # combine multiple cohorts
  sces <- lapply(sces, function(sce) {
    reducedDims(sce) <- reducedDims(sce)[reducedDimNames(sce) %in% dimred_to_keep]
    sce
  })
  sce <- Reduce(cbind, sces)
  rownames(sce) <- rowData(sce)[['Symbol']]
  
  # redo PCA, Harmony, UMAP, and TSNE; correcting for cohort level batch effect
  if (length(sces) > 1 & ncol(sce) > 20) {sce <- do_dimred(sce, harmonize = T, batch_col = "cohort", n_cores = snakemake@threads)}
  
  saveRDS(sces, snakemake@output[[paste0('scelist_', group)]])
  saveRDS(sce, snakemake@output[[paste0('sce_', group)]])
  write_tsv(sample_cell_counts, snakemake@output[[paste0('cells_per_sample_', group)]])

})


