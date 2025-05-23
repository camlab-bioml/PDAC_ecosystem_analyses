# Load libraries
suppressPackageStartupMessages({
        library(magrittr)
        library(tidyverse)
        library(sjstats)
        library(tidyr)
        library(stringr)
        library(SingleCellExperiment)
        library(scuttle)
        library(reticulate)
        library(singleCellTK)
        library(fastcluster)
        library(BiocParallel)
        library(GSVA)
        library(survival)
        # library(survminer)
        library(readxl)
})

# load TCGA expression data
PAAD.mRNA <- readRDS(snakemake@input[["paad_mrna_data"]])

## remove non PDAC samples
sample_correct <- read_xlsx(snakemake@input[["paad_sample_list"]], skip = 1)

keep_sample <- sample_correct[sample_correct[, 3] == "PDAC", 1][[1]]
keep_sample <- keep_sample[!is.na(keep_sample)]

samples_to_keep <- unlist(sapply(keep_sample, function(k) grep(k, PAAD.mRNA$gdc_cases.samples.submitter_id, value = T)))

PAAD.mRNA <- PAAD.mRNA[, PAAD.mRNA$gdc_cases.samples.submitter_id %in% samples_to_keep]
table(PAAD.mRNA$gdc_cases.samples.sample_type)
# table(PAAD.mRNA$gdc_cases.samples.submitter_id)
# colnames(PAAD.mRNA) <- PAAD.mRNA$gdc_cases.samples.submitter_id
head(colnames(PAAD.mRNA))

rm(keep_sample, samples_to_keep)

## normalize/transform expression
PAAD.dds <- DESeqDataSetFromMatrix(
        countData = PAAD.mRNA@assays$data$counts,
        colData = colData(PAAD.mRNA),
        design = ~gdc_cases.samples.submitter_id
)
PAAD.dds <- vst(PAAD.dds)

PAAD.mRNA.vst <- PAAD.dds@assays@data@listData[[1]]
colnames(PAAD.mRNA.vst) <- PAAD.mRNA$gdc_cases.samples.submitter_id

rm(PAAD.dds)

## update gene names from ENSEMBL ID to gene symbols
grch38 <- annotables::grch38

rownames(PAAD.mRNA.vst) <- gsub("\\.[1-9]|\\.[1-9][1-9]", "", rownames(PAAD.mRNA.vst))
PAAD.mRNA.vst <- PAAD.mRNA.vst[intersect(rownames(PAAD.mRNA.vst), grch38$ensgene), ]

rownames(PAAD.mRNA.vst) <- plyr::mapvalues(rownames(PAAD.mRNA.vst),
        from = grch38$ensgene,
        to = grch38$symbol,
        warn_missing = F
)

head(rownames(PAAD.mRNA.vst))
print(paste("Number of genes in PAAD mRNA data:", nrow(PAAD.mRNA.vst)))
print("Number of genes in PAAD mRNA data with missing gene symbols:")
table(is.na(rownames(PAAD.mRNA.vst)))

## transpose exprs matrix
PAAD.mRNA.vst <- PAAD.mRNA.vst |>
        t() |>
        as.data.frame()

## remove duplicated genes
print("Number of all genes:")
length(names(PAAD.mRNA.vst))
print("Number of unique genes:")
length(unique(names(PAAD.mRNA.vst)))

PAAD.mRNA.vst <- subset(PAAD.mRNA.vst, select=which(!duplicated(names(PAAD.mRNA.vst))))

## remove non expressed genes
print("Number of genes with uniform expression:")
table(apply(PAAD.mRNA.vst, 2, function(a) length(unique(a)) == 1))

PAAD.mRNA.vst <- PAAD.mRNA.vst[, !apply(PAAD.mRNA.vst, 2, function(a) length(unique(a)) == 1)]

## rownames (sample IDs) to column
PAAD.mRNA.vst <- PAAD.mRNA.vst |>
        rownames_to_column("bcr_patient_barcode")

saveRDS(PAAD.mRNA.vst, file = snakemake@output[["tidy_paad_mrna_data"]])