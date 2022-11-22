suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(rliger)
  library(tidyr)
  library(readr)
})

sces <- readRDS(snakemake@input[['scelist']])
var.genes <- read_tsv(snakemake@input[['var_genes_common']])

print("number of cells in each cohort: ")
print(sapply(sces, ncol))

# set optional knn_k for quantile_norm to the number of cells in the smallest cohort
min.cells <- sapply(sces, ncol) %>% min()
knn_k <- min.cells - 1

k = snakemake@params[['k']] %>% as.numeric()
lambda = snakemake@params[['Lambda']] %>% as.numeric()
seed = snakemake@params[['seed']] %>% as.numeric()

print("parameters for this LIGER run:")
print(paste0("seed: ", seed))
print(paste0("lambda: ", lambda))
print(paste0("k: ", k))

sces <- lapply(sces, function(sce) {
  sce[!grepl("^RP[LS]|^MT-", rownames(sce)),]
})

count_mtxlist <- lapply(sces, counts)

# Liger remove non-expressing genes in each cohort
liger <- createLiger(count_mtxlist,
                     take.gene.union = F,
                     remove.missing = T)

liger <- normalize(liger)

# use the common variable gene list between discovery and validation for signature extraction
# liger <- selectGenes(liger, 
#                      var.thresh = 0.1,
#                      alpha.thresh = 0.99,
#                      combine = "union",
#                      unshared = F, 
#                      unshared.datasets = NULL, 
#                      unshared.thresh= NULL,
#                      do.plot = F)
liger@var.genes <- var.genes$genes

# Liger removes cells with no gene expression
liger <- scaleNotCenter(liger,
                        remove.missing = T)

liger <- optimizeALS(liger, 
                     k = k, 
                     lambda = lambda, 
                     rand.seed = seed)
liger <- quantile_norm(liger, ref_dataset = NULL, knn_k = min(knn_k, 20))
liger <- louvainCluster(liger)
liger <- runUMAP(liger, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)

saveRDS(liger, file = snakemake@output[['liger']])


