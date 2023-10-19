suppressPackageStartupMessages({
        library(magrittr)
        library(tidyverse)
        library(sjstats)
        library(stringr)
        library(here)
        library(SingleCellExperiment)
        library(scater)
        library(scran)
})

sce.tumor <- readRDS(snakemake@input[['sce_tumor']])
sce.normal <- readRDS(snakemake@input[['sce_normal']])

table(sce.tumor$sample)
sce.tumor.subsampled <- sce.tumor
# sce.tumor.subsampled <- sce.tumor[, grepl("CRR241798|CRR034499|CRR034500|CRR034501|CRR034503|CRR034506|CRR034507|CRR034510|CRR034511|CRR034516|CRR034519", sce.tumor$sample)]
sce.normal <- sce.normal[,sample(colnames(sce.normal), 5000)]
table(colData(sce.tumor.subsampled)[[snakemake@params[["celltype_label_field"]]]])

sce.normal$cell_type <- paste0("Normal_", snakemake@wildcards[["celltype"]])
sce.tumor.subsampled$cell_type <- colData(sce.tumor.subsampled)[[snakemake@params[["celltype_label_field"]]]]

table(sce.tumor.subsampled$cell_type)
sce.tumor.subsampled$cell_type <- paste0("Tumour_", snakemake@wildcards[["celltype"]])
table(sce.tumor.subsampled$cell_type)

common.coldata <- intersect(names(colData(sce.normal)), names(colData(sce.tumor.subsampled)))

colData(sce.normal) <- colData(sce.normal)[, common.coldata]
colData(sce.tumor.subsampled) <- colData(sce.tumor.subsampled)[, common.coldata]

reducedDims(sce.normal) <- NULL
reducedDims(sce.tumor.subsampled) <- NULL

head(rownames(sce.tumor.subsampled))
head(rownames(sce.normal))

rownames(sce.normal) <- paste(rowData(sce.normal)[["Symbol"]], rowData(sce.normal)[["ensembl_id"]], sep = "_")
rownames(sce.tumor.subsampled) <- paste(rowData(sce.tumor)[["Symbol"]], rowData(sce.tumor.subsampled)[["ensembl_id"]], sep = "_")

common.gene <- intersect(rownames(sce.normal), rownames(sce.tumor.subsampled))
sce.normal <- sce.normal[common.gene, ]
sce.tumor.subsampled <- sce.tumor.subsampled[common.gene, ]

rownames(sce.normal) <- str_split(rownames(sce.normal), pattern = "_", simplify = T)[, 1]
rownames(sce.tumor.subsampled) <- str_split(rownames(sce.tumor.subsampled), pattern = "_", simplify = T)[, 1]


rowData(sce.normal)[["seurat_variableFeatures_vst_varianceStandardized"]] <- NULL
rowData(sce.normal)[["seurat_variableFeatures_vst_mean"]] <- NULL
rowData(sce.normal)[["chr"]] <- NULL
rowData(sce.normal)[["gene_start"]] <- NULL
rowData(sce.normal)[["gene_end"]] <- NULL
rowData(sce.normal)[["gene_strand"]] <- NULL

rowData(sce.tumor.subsampled)[["seurat_variableFeatures_vst_varianceStandardized"]] <- NULL
rowData(sce.tumor.subsampled)[["seurat_variableFeatures_vst_mean"]] <- NULL
rowData(sce.tumor.subsampled)[["chr"]] <- NULL
rowData(sce.tumor.subsampled)[["gene_start"]] <- NULL
rowData(sce.tumor.subsampled)[["gene_end"]] <- NULL
rowData(sce.tumor.subsampled)[["gene_strand"]] <- NULL

sce.combined <- cbind(sce.tumor.subsampled, sce.normal)

saveRDS(sce.combined, file = snakemake@output[['sce_combined']])

write_csv(data.frame(sample = unique(sce.tumor.subsampled$sample)), file = snakemake@output[['sample_list']])
