suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(rliger)
  library(sjstats)
  library(tidyr)
  library(stringr)
})

discovery <- readRDS(snakemake@input[['liger_dis']])
validation <- readRDS(snakemake@input[['liger_val']])
sce_dis <- readRDS(snakemake@input[['sce_dis']])
sce_val <- readRDS(snakemake@input[['sce_val']])

# get signature loading matrix
h.discovery <- discovery@H.norm
h.validation <- validation@H.norm

colnames(h.discovery) <- paste0("discovery ", seq(ncol(h.discovery)))
colnames(h.validation) <- paste0("validation ", seq(ncol(h.validation)))

h.discovery <- h.discovery %>% as.data.frame()
h.validation <- h.validation %>% as.data.frame()

h.discovery$cohort <- sce_dis$cohort
h.discovery$sample <- sce_dis$sample
h.discovery$cell_id <- rownames(h.discovery)

h.validation$cohort <- sce_val$cohort
h.validation$sample <- sce_val$sample
h.validation$cell_id <- rownames(h.validation)

# get gene loading matrix 
w.discovery <- discovery@W
w.validation <- validation@W

common.genes <- intersect(colnames(w.discovery), colnames(w.validation))
print(paste0("number of genes in discovery signatures for ", snakemake@wildcards[['subtype']], ": ", length(colnames(w.discovery))))
print(paste0("number of genes in validation signatures for ", snakemake@wildcards[['subtype']], ": ", length(colnames(w.validation))))
print(paste0("number of common genes between discovery signatures and validation signatures for ", snakemake@wildcards[['subtype']], ": ", length(common.genes)))

rownames(w.discovery) <- paste0("discovery ", seq(nrow(w.discovery)))
rownames(w.validation) <- paste0("validation ", seq(nrow(w.validation)))

w.discovery <- w.discovery %>% t() %>% as.data.frame()
w.validation <- w.validation %>% t() %>% as.data.frame()

w.discovery$gene <- rownames(w.discovery)
w.validation$gene <- rownames(w.validation)

# union join discovery and validation group results
w <- full_join(w.discovery, w.validation, by = "gene")
rownames(w) <- w$gene

# save the loading matrices
write_tsv(w, file = snakemake@output[['gene_loading_mtx']])
write_tsv(h.discovery, file = snakemake@output[['sig_loading_mtx_dis']])
write_tsv(h.validation, file = snakemake@output[['sig_loading_mtx_val']])

