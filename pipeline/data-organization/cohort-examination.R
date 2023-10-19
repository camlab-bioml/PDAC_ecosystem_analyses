suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(scater)
  library(scran)
  library(singleCellTK)
  library(BiocParallel)
  library(sjstats)
})

set.seed(123L)

celltypes <- snakemake@params[['celltypes']]

# load sces
sce.list <- lapply(snakemake@input[['sces']], function(sce.path) {
  readRDS(sce.path)
})
names(sce.list) <- celltypes

print("SCEs successfully loaded")

lapply(sce.list, function(sce) {
  head(rownames(sce))
  dim(sce)
})

## Some sanity checks
stopifnot(all(sapply(sce.list, function(sce) {length(unique(sce$cohort))}) > 1))

# group cohorts
groups <- list(
  discovery = snakemake@params[["cohorts_discovery"]],
  validation = snakemake@params[["cohorts_validation"]]
)

cohorts <- groups[[snakemake@wildcards[["group"]]]]

# subset to common genes across cohorts
sce.list <- lapply(sce.list, function(sce) {
  rownames(sce) <- paste(rowData(sce)$ensembl_id, rownames(sce), sep = "_")
  sce
})
genes.common <- Reduce(intersect, lapply(sce.list, rownames))
sce.list <- lapply(sce.list, function(sce) sce[genes.common, ])

# cbind celltype sces
reddims.common <- Reduce(intersect, lapply(sce.list, reducedDimNames))
coldata.list <- lapply(sce.list, colData)
coldata.common <- Reduce(intersect, lapply(coldata.list, names))

## subset to commone metadata
sce.list <- lapply(sce.list, function(sce) {
  colData(sce) <- colData(sce) %>%
    as.data.frame() %>%
    select(all_of(coldata.common)) %>%
    DataFrame()
  reducedDims(sce) <- reducedDims(sce)[reddims.common]
  reducedDims(sce)["HARMONY"] <- NULL
  reducedDims(sce)["PCA_ON_HARMONY"] <- NULL
  reducedDims(sce)["PCA_on_Harmony"] <- NULL
  reducedDims(sce)["Harmony"] <- NULL

  rowData(sce)[["seurat_variableFeatures_vst_varianceStandardized"]] <- NULL
  rowData(sce)[["seurat_variableFeatures_vst_mean"]] <- NULL
  rowData(sce)[["chr"]] <- NULL
  rowData(sce)[["gene_start"]] <- NULL
  rowData(sce)[["gene_end"]] <- NULL
  rowData(sce)[["gene_strand"]] <- NULL

  sce
})

sce <- Reduce(cbind, sce.list)
sce$celltype <- sce$singler.pruned.label

print(paste0("Celltypes successfully combined into one SCE for ", snakemake@wildcards[["group"]]))

# subset combined sce to split cohorts
sce.list <- lapply(cohorts, function(cohort) {
  sce <- sce[, colData(sce)[["cohort"]] == cohort]
  sce
})

print(paste0("Cohorts successfully split for ", snakemake@wildcards[["group"]]))

# compute Seurat normalized counts
sce.list <- lapply(sce.list, function(sce) {
  sce <- SingleCellExperiment(assays = assays(sce), rowData = rowData(sce), colData = colData(sce), reducedDims = reducedDims(sce), altExps = NULL)
  sce <- scuttle::logNormCounts(sce, BPPARAM = MulticoreParam(snakemake@threads))
  sce <- runSeuratNormalizeData(sce, useAssay = "counts")
  sce
})

print(paste0("Seurat normalized counts and logcounts successfully computed for ", snakemake@wildcards[["group"]]))

# sample sces
sce.list <- lapply(sce.list, function(sce) {
  sce <- sce[, sample(ncol(sce), 1000)]
  sce
})

print(paste0("SCEs successfully sampled for ", snakemake@wildcards[["group"]]))

# subset to common genes across celltypes
genes.common <- Reduce(intersect, lapply(sce.list, rownames))
sce.list <- lapply(sce.list, function(sce) sce[genes.common, ])

# cbind celltype sces
reddims.common <- Reduce(intersect, lapply(sce.list, reducedDimNames))
coldata.list <- lapply(sce.list, colData)
coldata.common <- Reduce(intersect, lapply(coldata.list, names))

## subset to commone metadata
sce.list <- lapply(sce.list, function(sce) {
  colData(sce) <- colData(sce) %>%
    as.data.frame() %>%
    select(all_of(coldata.common)) %>%
    DataFrame()
  reducedDims(sce) <- reducedDims(sce)[reddims.common]
  reducedDims(sce)["HARMONY"] <- NULL
  reducedDims(sce)["PCA_ON_HARMONY"] <- NULL
  reducedDims(sce)["PCA_on_Harmony"] <- NULL
  reducedDims(sce)["Harmony"] <- NULL

  rowData(sce)[["seurat_variableFeatures_vst_varianceStandardized"]] <- NULL
  rowData(sce)[["seurat_variableFeatures_vst_mean"]] <- NULL
  rowData(sce)[["chr"]] <- NULL
  rowData(sce)[["gene_start"]] <- NULL
  rowData(sce)[["gene_end"]] <- NULL
  rowData(sce)[["gene_strand"]] <- NULL

  sce
})

sce <- Reduce(cbind, sce.list)
sce$celltype <- sce$singler.pruned.label

print(paste0("Cohorts successfully combined into one SCE for ", snakemake@wildcards[["group"]]))


# save sce
saveRDS(sce, file = snakemake@output[['sce']])
print("SCE successfully saved")

