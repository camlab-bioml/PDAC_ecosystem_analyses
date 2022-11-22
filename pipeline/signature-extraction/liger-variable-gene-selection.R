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

sces_dis <- readRDS(snakemake@input[['sces_dis']])
sces_val <- readRDS(snakemake@input[['sces_val']])

exprs.var.thres <- snakemake@params[['exprs_var_thres']] %>% as.numeric()

print("number of cells in each discovery cohort: ")
print(sapply(sces_dis, ncol))

print("number of cells in each validation cohort: ")
print(sapply(sces_val, ncol))

# remove mito and ribo genes from the sces
sces_dis <- lapply(sces_dis, function(sce) {
  sce[!grepl("^RP[LS]|^MT-", rownames(sce)),]
})

sces_val <- lapply(sces_val, function(sce) {
  sce[!grepl("^RP[LS]|^MT-", rownames(sce)),]
})

count_mtxlist_dis <- lapply(sces_dis, counts)
count_mtxlist_val <- lapply(sces_val, counts)

# LIGER removes non-expressing genes in each cohort
liger.dis <- createLiger(count_mtxlist_dis,
                         take.gene.union = F,
                         remove.missing = T)

liger.val <- createLiger(count_mtxlist_val,
                         take.gene.union = F,
                         remove.missing = T)

liger.dis <- normalize(liger.dis)
liger.val <- normalize(liger.val)

# find variable genes in the discovery group, use this var.genes list for signature extraction in the next step
liger.dis <- selectGenes(liger.dis, 
                         var.thresh = exprs.var.thres,
                         alpha.thresh = 0.99,
                         combine = "union",
                         unshared = F, 
                         unshared.datasets = NULL, 
                         unshared.thresh= NULL,
                         do.plot = F)

liger.val <- selectGenes(liger.val, 
                         var.thresh = exprs.var.thres,
                         alpha.thresh = 0.99,
                         combine = "union",
                         unshared = F, 
                         unshared.datasets = NULL, 
                         unshared.thresh= NULL,
                         do.plot = F)

# get common variable genes between discovery group and validation group 
genes.val.common <- lapply(liger.val@norm.data, function(mtx) {
  mtx@Dimnames[[1]]
})
if(length(genes.val.common) == 1) {
  genes.val.common <- genes.val.common[[1]]
} else {
  genes.val.common <- Reduce(intersect, genes.val.common)
}

print("number of common genes in the validation cohorts: ")
print(genes.val.common %>% length())

var.genes.common <- intersect(genes.val.common, liger.dis@var.genes)

print("number of common variable genes between discovery and validation: ")
print(var.genes.common %>% length())

# save the gene lists
write_tsv(data_frame(genes = liger.dis@var.genes), file = snakemake@output[['var_genes_dis']])
write_tsv(data_frame(genes = liger.val@var.genes), file = snakemake@output[['var_genes_val']])
write_tsv(data_frame(genes = var.genes.common), file = snakemake@output[['var_genes_common']])


