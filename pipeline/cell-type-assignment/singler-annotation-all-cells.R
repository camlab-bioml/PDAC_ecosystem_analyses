suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(scuttle)
  library(SingleR)
  library(BiocParallel)
})

source("utils/myhelpers.R")

# load data
sceQuery <- readRDS(snakemake@input[['sce_query']])
sceRef <- readRDS(snakemake@input[['sce_ref']])
sceRef <- Seurat::as.SingleCellExperiment(sceRef, assay = "RNA")
print("SCEs loaded")

# change rownames of the count matrices of the sce to remove ENSEMBL IDs so that the rownames matches the SingleR reference rownames
print("changing rownames of the query SCE from: ")
head(rownames(sceQuery))
print("to match the reference SCE rownames: ")
head(rownames(sceRef))
rowData(sceQuery)[['ensembl_id']] <- str_split(rownames(sceQuery), "_", simplify = T)[,1]

gene_names <- rowData(sceQuery)[["Symbol"]] %>% make.unique()
rownames(sceQuery) <- gene_names

# SingleR uses log-transformed expression from the reference dataset
stopifnot("logcounts" %in% assayNames(sceRef))

ref_label_field = snakemake@params[['ref_label_field']]
print("reference cell type label field:")
print(ref_label_field)
print("reference cell types:")
table(colData(sceRef)[[ref_label_field]])

# subset to cells from a specific project (Peng) in the reference
sceRef.to.use <- sceRef[,sceRef$Project == "CA001063"]

# rowData(sceRef.to.use)[["symbol"]] <- rownames(sceRef.to.use)
# rowData(sceRef.to.use)[['ensembl_id']] <- plyr::mapvalues(rowData(sceRef.to.use)[["symbol"]], 
#                                                           from = annotables::grch38$symbol, 
#                                                           to = annotables::grch38$ensgene,
#                                                           warn_missing = F)
# 
# rownames(sceRef.to.use) <- paste(rowData(sceRef.to.use)[["ensembl_id"]], rowData(sceRef.to.use)[["symbol"]], sep = "_")


print(paste0("SingleR running for ", snakemake@wildcards[['cohort']]))
common_genes = intersect(rownames(sceQuery), rownames(sceRef.to.use))
num_common_genes = common_genes %>% length()
print(paste0("number of genes in common between query and reference: ", num_common_genes))

# we want 5 markers per cell type on average
stopifnot(num_common_genes >= 5 * length(unique(colData(sceRef.to.use)[[ref_label_field]])))

# we want to remove mitochondrial and ribosomal genes from marker selection for SingleR
print("mitochondrial and ribosomal genes present in the common gene list:")
print(grep("^RP[LS]|^MT-", common_genes, value = T))
print("they will be removed from marker selection")

# many parameters should be tuned when running SingleR
pred <- SingleR(test = sceQuery,
                ref = sceRef.to.use,
                labels = colData(sceRef.to.use)[[ref_label_field]],
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

print(paste0("SingleR finished for ", snakemake@wildcards[["cohort"]]))

sceQuery$singler.first.label <- pred$first.labels
sceQuery$singler.pruned.label <- pred$pruned.labels
sceQuery$singler.label <- pred$labels
sceQuery$singler.best.score <- pred$tuning.scores$first

saveRDS(sceQuery, file = snakemake@output[["sce"]])
print(paste0("Updated SCE saved for ", snakemake@wildcards[["cohort"]]))
saveRDS(pred, file = snakemake@output[["dframe"]])





