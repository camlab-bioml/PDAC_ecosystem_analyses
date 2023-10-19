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

# get cell meta data
coldata.list <- lapply(sce.list, colData)
coldata.common <- Reduce(intersect, lapply(coldata.list, names))
coldata.list <- lapply(coldata.list, function(coldata) {
  coldata %>% as.data.frame() %>% select(all_of(coldata.common))
})
coldata <- Reduce(rbind, coldata.list) %>% as.data.frame()

# group cohorts
groups <- list(discovery = snakemake@params[['cohorts_discovery']],
               validation = snakemake@params[['cohorts_validation']])

cohorts <- groups[[snakemake@wildcards[["group"]]]]

# get metadata of interest
coldata <- coldata %>% 
  group_by(cohort, sample) %>%
  summarise(num_donor = n_distinct(sample),
            num_cell = n(),
            num_genes = mean(detected),
            num_UMI = mean(sum),
            percent_mito = mean(subsets_mito_percent),
            percent_ribo = mean(subsets_ribo_percent)) %>%
  ungroup() %>%
  group_by(cohort) %>%
  mutate(num_donor = n_distinct(sample)) %>%
  ungroup()

print("Metadata successfully computed")

# save metadata
saveRDS(coldata, file = snakemake@output[['metadata']])

# subset to common genes across celltypes
sce.list <- lapply(sce.list, function(sce) {
  rownames(sce) <- paste(rowData(sce)$ensembl_id, rownames(sce), sep = "_")
  sce
})
genes.common <- Reduce(intersect, lapply(sce.list, rownames))
sce.list <- lapply(sce.list, function(sce) sce[genes.common,])

# cbind celltype sces
reddims.common <- Reduce(intersect, lapply(sce.list, reducedDimNames))

## subset to commone metadata
sce.list <- lapply(sce.list, function(sce) {
  colData(sce) <- colData(sce) %>% as.data.frame() %>% select(all_of(coldata.common)) %>% DataFrame()
  reducedDims(sce) <- reducedDims(sce)[reddims.common]
  reducedDims(sce)["HARMONY"] <- NULL
  reducedDims(sce)["PCA_ON_HARMONY"] <- NULL
  reducedDims(sce)["PCA_on_Harmony"] <- NULL
  reducedDims(sce)["Harmony"] <- NULL
  #logcounts(sce) <- NULL
  
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

print(paste0("Celltypes successfully combined into one SCE for ", snakemake@wildcards[['group']]))

# run Seurat PCA and UMAP for better plots
set.seed(42)
sce <- SingleCellExperiment(assays = assays(sce), rowData = rowData(sce), colData = colData(sce), reducedDims = NULL, altExps = NULL)

##### remove ribo/mito genes #####
sce <- sce[!grepl("_(RP|MT)-", rownames(sce)),]

sce <- runSeuratNormalizeData(sce, useAssay = "counts")
sce <- runSeuratFindHVG(sce, useAssay = "counts")
#sce <- runSeuratIntegration(sce, batch = "cohort", useAssay = "counts", newAssayName = "SeuratIntegratedAssay", kAnchor = 5, kFilter = 200, kWeight = 100)
assay(sce, "counts") <- NULL
sce <- runSeuratScaleData(sce,
                          useAssay = "seuratNormData"
                          #useAssay = "SeuratIntegratedAssay"
                          )
#assay(sce, "SeuratNormData") <- NULL
#assay(sce, "SeuratIntegratedAssay") <- NULL
print("Available assays:")
assays(sce)
sce <- runSeuratPCA(sce, useAssay = "seuratNormData")
#sce <- runSeuratPCA(sce, useAssay = "seuratScaledData")
#sce <- runHarmony(sce, useReducedDim = "seuratPCA", batch = "cohort", reducedDimName = "HARMONY", theta = 2, sigma = 0.1, nComponents = 50)
sce <- runSeuratFindClusters(sce, useAssay = "seuratNormData", useReduction = "pca", dims = 10, algorithm = "louvain", resolution = 0.8)
print("Available reducedDims:")
reducedDims(sce)
ncol(sce)
length(unique(colnames(sce)))
#dim(reducedDim(sce, "HARMONY"))
#sce <- runSeuratUMAP(sce, nNeighbors = 30, externalReduction = Seurat::CreateDimReducObject(embeddings = reducedDim(sce, "HARMONY"), key = "HARMONY"), seed = 42, verbose = TRUE)
sce <- runSeuratUMAP(sce, nNeighbors = 30, useReduction = "pca", seed = 42, verbose = TRUE)
#sce <- runUMAP(sce, useReducedDim = "HARMONY") # leave everything else as defaults
#sce <- runUMAP(sce, useReducedDim = "seuratPCA") # leave everything else as defaults
sce



print(paste0("Seurat UMAP successfully computed for ", snakemake@wildcards[['group']]))
assayNames(sce)

# save sce
saveRDS(sce, file = snakemake@output[['sce']])

# get desired metadata
dim.red.plot = snakemake@params[['dim_red_plot']]

df.redim <- data.frame(Cell_ID = colnames(sce),
                       Sample = sce$sample,
                       Cohort = sce$cohort,
                       Cell_type = sce$celltype,
                       UMAP_1 = reducedDim(sce, dim.red.plot)[,1],
                       UMAP_2 = reducedDim(sce, dim.red.plot)[,2],
                       Seurat_cluster = sce$Seurat_louvain_Resolution0.8)

# tidyup metadata
print("Cell types:")
print(table(df.redim$Cell_type))

df.redim <- df.redim %>%
    # add column n with counts for each celltype
    add_count(Cell_type) %>% 
    # combine the cell type and count n into one column
    mutate(Cell_type = paste0(Cell_type, ' (', n, ')')) %>%
    # add column n with counts for each cohort
    add_count(Cohort) %>% 
    # combine the cohort and count n into one column
    mutate(Cohort = paste0(Cohort, ' (', nn, ')')) %>%
    # set 'NA (*)' labels to NA
    mutate(Cell_type = ifelse(grepl("NA", Cell_type), NA, Cell_type))

print(paste0("Metadata for ", dim.red.plot, " successfully summarized"))

print(table(df.redim$Cell_type) %>% sort(decreasing = T))
print(table(df.redim$Cohort) %>% sort(decreasing = T))

# adding number of cells per cohort in the legend
cohorts.with.num.cell <- df.redim$Cohort %>% unique()

cohorts.with.num.cell
factor(str_split(cohorts.with.num.cell, " ", simplify = T)[,1], levels = cohorts) %>% rank()
cohorts.with.num.cell[order(factor(str_split(cohorts.with.num.cell, " ", simplify = T)[, 1], levels = cohorts) %>% rank())]
cohorts.with.num.cell

# save df.redim.list
saveRDS(df.redim, file = snakemake@output[['dimred']])

print(paste0("UMAP coordinates and cell metadata successfully saved for ", snakemake@wildcards[['group']]))


