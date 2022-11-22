suppressPackageStartupMessages({
  library(magrittr)
  library(SingleCellExperiment)
  library(scuttle)
  library(scater)
  library(SingleR)
})

source("utils/myhelpers.R")

sceQuery <- readRDS(snakemake@input[['sce']])
sceRef <- readRDS(snakemake@input[['ref']])

# SingleR uses log-transformed expression from the reference dataset
stopifnot("logcounts" %in% assayNames(sceRef))

ref_label_field = snakemake@params[['ref_label_field']]
print("reference cell type label field:")
print(ref_label_field)

common_genes = intersect(rownames(sceQuery), rownames(sceRef))
num_common_genes = common_genes %>% length()
print(paste0("number of genes in common between query and reference: ", num_common_genes))

# we want 5 markers per cell type on average
stopifnot(num_common_genes >= 5 * length(unique(colData(sceRef)[[ref_label_field]])))

# we want to remove mitochondrial and ribosomal genes from marker selection for SingleR
print("mitochondrial and ribosomal genes present in the common gene list:")
print(grep("^RP[LS]|^MT-", common_genes, value = T))
print("they will be removed from marker selection")

# many parameters should be tuned when running SingleR 
pred <- SingleR(test = sceQuery, 
                ref = sceRef, 
                labels = colData(sceRef)[[ref_label_field]], 
                restrict = filter_genes(common_genes, "^RP[LS]|^MT-", return.gnames = T),
                genes = "de",
                sd.thresh = 1,
                de.method = "classic",
                de.n = NULL, 
                de.args = list(), 
                aggr.ref = T,
                aggr.args = list(power = 0.5),
                recompute = T,
                quantile = 0.8,
                fine.tune = T,
                tune.thresh = 0.05,
                prune = T,
                assay.type.test = "logcounts",
                assay.type.ref = "logcounts",
                BPPARAM = BiocParallel::MulticoreParam(as.numeric(snakemake@threads)))

sceQuery$singler.first.label <- pred$first.labels
sceQuery$singler.pruned.label <- pred$pruned.labels
sceQuery$singler.label <- pred$labels
sceQuery$singler.best.score <- pred$tuning.scores$first

saveRDS(sceQuery, file = snakemake@output[['sce']])
saveRDS(pred, file = snakemake@output[['dframe']])

scores.no.fine.tune <- pred$scores
rownames(scores.no.fine.tune) <- colnames(sceQuery)
readr::write_tsv(as.data.frame(scores.no.fine.tune), file = snakemake@output[['tsv']])

