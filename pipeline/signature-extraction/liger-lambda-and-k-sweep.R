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

scelist <- readRDS(snakemake@input[['scelist']])
scelist <- lapply(scelist, function(sce) {
  sce[!grepl("^RP[LS]|^MT-", rownames(sce)),]
})
count.mtx.list <- lapply(scelist, counts)

print("number of cells in each cohort: ")
print(sapply(scelist, ncol))

# set optional knn_k for quantile_norm to the number of cells in the smallest cohort
min.cells <- sapply(scelist, ncol) %>% min()
knn_k <- min.cells - 1

seed <- snakemake@params[['seed']] %>% as.numeric()
lambda <- snakemake@params[['Lambda']] %>% as.numeric()
k <- snakemake@params[['k']] %>% as.numeric()

if (min.cells < k) {
  print(paste0("selected k: ", k))
  print(paste0("smallest cohort has lower number of cells than selected k, setting k and knn_k to: ", min.cells))
  k <- min.cells
}

print("parameters for this LIGER run:")
print(paste0("seed: ", seed))
print(paste0("lambda: ", lambda))
print(paste0("k: ", k))

# Liger remove non-expressing genes in each cohort
liger <- createLiger(count.mtx.list, take.gene.union = F, remove.missing = T) %>% normalize()
cohorts <- names(liger@raw.data)

rm(scelist)
rm(count.mtx.list)

liger <- selectGenes(liger, 
                     var.thresh = 0.1,
                     alpha.thresh = 0.99,
                     combine = "union",
                     unshared = F, 
                     unshared.datasets = NULL, 
                     unshared.thresh= NULL,
                     do.plot = F)

# Liger removes cells with no gene expression
liger <- scaleNotCenter(liger, remove.missing = T)

# run Liger with a set of parameters
liger <- optimizeALS(liger, lambda = lambda, k = k, rand.seed = seed, use.unshared = F, max.iters = 30, thresh = 1e-10)
liger <- quantile_norm(liger, ref_dataset = NULL, knn_k = min(knn_k, 20))
liger <- louvainCluster(liger)

# calculate alignment score and agreement score
alignment = calcAlignment(liger, rand.seed = seed, by.dataset = T)
alignment_all = calcAlignment(liger, rand.seed = seed, by.dataset = F)
alignment_df <- c(alignment, alignment_all) %>% as.matrix() %>% t() %>% as.data.frame()
names(alignment_df) <- c(cohorts, "overall")
rownames(alignment_df) <- c("alignment")
alignment_df$parameter <- paste0("Seed:", seed, "_Lambda:", lambda, "_K:",k)

agreement = calcAgreement(liger, dr.method = "NMF", ndims = k, k = min(knn_k, 15), use.aligned = T, rand.seed = seed, by.dataset = T)
agreement_all = calcAgreement(liger, dr.method = "NMF", ndims = k, k = min(knn_k, 15), use.aligned = T, rand.seed = seed, by.dataset = F)
agreement_df <- c(agreement, agreement_all) %>% as.matrix() %>% t() %>% as.data.frame()
names(agreement_df) <- c(cohorts, "overall")
rownames(agreement_df) <- c("agreement")
agreement_df$parameter <- paste0("Seed:", seed, "_Lambda:", lambda, "_K:",k)

metrics_df <- bind_rows(alignment_df, agreement_df)
metrics_df$metric <- rownames(metrics_df)

# calculate contributions from each matrix
normloadings = calcNormLoadings(liger)
names(normloadings) <- c("factor", "W_loadings", paste0(cohorts, "_loadings"))

write_tsv(metrics_df, file = snakemake@output[['metrics']])
write_tsv(normloadings, file = snakemake@output[['factor_contributions']])




