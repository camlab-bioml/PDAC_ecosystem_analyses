suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(Azimuth)
  library(SingleCellExperiment)
  library(magrittr)
  library(stringr)
  library(tidyverse)
})


sce <- readRDS(snakemake@input[['sce']])

# change rownames of the count matrices of the sce to remove ENSEMBL IDs so that the rownames matches the Azimuth reference rownames
rowData(sce)[['ensembl_id']] <- str_split(rownames(sce), "_", simplify = T)[,1]

gene_names <- rowData(sce)[["Symbol"]] %>% make.unique()
rownames(sce) <- gene_names

exprs_mtx <- counts(sce)
rownames(exprs_mtx) <- gene_names
counts(sce) <- exprs_mtx

exprs_mtx <- logcounts(sce)
rownames(exprs_mtx) <- gene_names
logcounts(sce) <- exprs_mtx

rm(exprs_mtx)

seu <- as.Seurat(sce, counts = "counts", data = "logcounts", project = "SingleCellExperiment")
rm(sce)

# Run Azimuth annotation on Seurat object
# Availability of annotation levels (l1 or l2/l3) is different when using different references, check Azimuth website for details
seu <- RunAzimuth(query = seu, 
                  assay = "originalexp",
                  reference = snakemake@params[['azimuth_ref']],
                  annotation.levels = c(snakemake@params[['annot_level']]))

# Coerce back to sce object
sce_annot <- as.SingleCellExperiment(seu, assay = "originalexp")
rowData(sce_annot) <- seu@assays$originalexp@meta.features
colData(sce_annot) <- colData(sce_annot) %>% 
  as.data.frame() %>% 
  select(!contains(c("ident", "_originalexp", "_refAssay", "_RNA", "percent.mt"))) %>%
  as("DataFrame")

# pull out predicted cell type annotation scores for each cell
# 20220714: be careful about annotation level naming for different references, e.g. for pancreas reference there is annotation.l1, for pbmc reference there are celltype.l1 (or l2/l3)
score_mtx <- seu@assays[[paste0("prediction.score.", snakemake@params[['annot_level']])]]@data %>% t() %>% as.data.frame()
score_mtx$cell_id <- rownames(score_mtx)
score_mtx$cell_type <- colData(sce_annot)[[paste0("predicted.", snakemake@params[['annot_level']])]]

saveRDS(seu, file = snakemake@output[['seu']])
saveRDS(sce_annot, file = snakemake@output[['sce']])
write_tsv(score_mtx, file = snakemake@output[['assignment_score']])




