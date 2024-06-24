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
        library(DESeq2)
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

# load PanCuRx + COMPASS expression and clinical data
toronto.mRNA <- readRDS(snakemake@input[["toronto_bulk_data"]])
compass.clin <- read_tsv(snakemake@input[["tidy_compass_clin_data"]])
pancurx.clin.new <- read_tsv(snakemake@input[["tidy_pancurx_clin_data"]])

## subset to compass/pancurx samples
compass.mRNA <- toronto.mRNA[, toronto.mRNA$PCSI_ID %in% compass.clin$PCSI_ID]
pancurx.mRNA.new <- toronto.mRNA[, toronto.mRNA$sample_id %in% pancurx.clin.new$sample_id]

print("Number of samples in COMPASS bulk data:")
length(unique(compass.mRNA$sample_id))
print("Number of patients in COMPASS bulk data:")
length(unique(compass.mRNA$PCSI_ID))

print("Tissue status in COMPASS bulk data:")
is.na(compass.mRNA$tissue) |> table()
table(compass.mRNA$tissue)

print("Tissue status in PanCuRx bulk data:")
is.na(pancurx.mRNA.new$tissue) |> table()
table(pancurx.mRNA.new$tissue)

print("Number of Pa patients in COMPASS bulk data:")
length(unique(compass.mRNA[, compass.mRNA$tissue == "Pa"]$PCSI_ID))
print("Number of Lv patients in COMPASS bulk data:")
length(unique(compass.mRNA[, compass.mRNA$tissue == "Lv"]$PCSI_ID))

## subset to samples at pancreas/liver
pancurx.mRNA.logTPM <- pancurx.mRNA.new

compass.Pa.mRNA <- compass.mRNA[, compass.mRNA$tissue == "Pa"]
compass.Lv.mRNA <- compass.mRNA[, compass.mRNA$tissue == "Lv"]

compass.mRNA.logTPM <- compass.mRNA
compass.Pa.mRNA.logTPM <- compass.Pa.mRNA
compass.Lv.mRNA.logTPM <- compass.Lv.mRNA

## subset to genes without ENSEMBL ID
pancurx.mRNA.logTPM <- pancurx.mRNA.logTPM[!grepl("_ENSG", rownames(pancurx.mRNA.logTPM)), ]
compass.mRNA.logTPM <- compass.mRNA.logTPM[!grepl("_ENSG", rownames(compass.mRNA.logTPM)), ]
compass.Pa.mRNA.logTPM <- compass.Pa.mRNA.logTPM[!grepl("_ENSG", rownames(compass.Pa.mRNA.logTPM)), ]
compass.Lv.mRNA.logTPM <- compass.Lv.mRNA.logTPM[!grepl("_ENSG", rownames(compass.Lv.mRNA.logTPM)), ]

## transpose exprs matrix
pancurx.mRNA.logTPM <- assay(pancurx.mRNA.logTPM, "logTPM")
compass.mRNA.logTPM <- assay(compass.mRNA.logTPM, "logTPM")
compass.Pa.mRNA.logTPM <- assay(compass.Pa.mRNA.logTPM, "logTPM")
compass.Lv.mRNA.logTPM <- assay(compass.Lv.mRNA.logTPM, "logTPM")

pancurx.mRNA.logTPM <- pancurx.mRNA.logTPM |>
        t() |>
        as.data.frame() |>
        rownames_to_column("sample_id")
compass.mRNA.logTPM <- compass.mRNA.logTPM |>
        t() |>
        as.data.frame() |>
        rownames_to_column("sample_id")
compass.Pa.mRNA.logTPM <- compass.Pa.mRNA.logTPM |>
        t() |>
        as.data.frame() |>
        rownames_to_column("sample_id")
compass.Lv.mRNA.logTPM <- compass.Lv.mRNA.logTPM |>
        t() |>
        as.data.frame() |>
        rownames_to_column("sample_id")

## remove duplicated genes
pancurx.mRNA.logTPM <- subset(pancurx.mRNA.logTPM, select = which(!duplicated(names(pancurx.mRNA.logTPM))))
compass.mRNA.logTPM <- subset(compass.mRNA.logTPM, select=which(!duplicated(names(compass.mRNA.logTPM))))
compass.Pa.mRNA.logTPM <- subset(compass.Pa.mRNA.logTPM, select=which(!duplicated(names(compass.Pa.mRNA.logTPM))))
compass.Lv.mRNA.logTPM <- subset(compass.Lv.mRNA.logTPM, select=which(!duplicated(names(compass.Lv.mRNA.logTPM))))

saveRDS(pancurx.mRNA.logTPM, file = snakemake@output[["tidy_pancurx_mrna_data"]])
saveRDS(compass.mRNA.logTPM, file = snakemake@output[["tidy_compass_mrna_data"]])
saveRDS(compass.Pa.mRNA.logTPM, file = snakemake@output[["tidy_compass_pa_mrna_data"]])
saveRDS(compass.Lv.mRNA.logTPM, file = snakemake@output[["tidy_compass_lv_mrna_data"]])