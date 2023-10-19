suppressPackageStartupMessages({
    library(infercnv)
    library(SingleCellExperiment)
    library(tibble)
    library(dplyr)
    library(tidyr)
    library(argparse)
})

set.seed(123L)

sce <- readRDS(snakemake@input[["sce_combined"]])

sub_sample <- sce[, colData(sce)[["sample"]] == snakemake@wildcards[["sample"]] ]

## Some sanity checks
stopifnot(ncol(sub_sample) > 10)

sce_tumour <- sub_sample[, colData(sub_sample)[["cell_type"]] == paste0("Tumour_", snakemake@wildcards[["celltype"]])]

sce_normal <- sce[, colData(sce)[["cell_type"]] == paste0("Normal_", snakemake@wildcards[["celltype"]])]

if(ncol(sce_normal) > 3000) {
    sce_normal <- sce_normal[, sample(ncol(sce_normal), 3000)]
}

patient_sce <- cbind(sce_tumour, sce_normal)

print(dim(sce_tumour))
print(dim(sce_normal))
print(dim(patient_sce))

cts <- assay(patient_sce, 'counts')

gene_symbols <- gsub("ENSG[0-9]*-", "", rownames(sce))
count_mat <- as.matrix(assay(patient_sce, 'counts'))
cts <- rowsum(count_mat, gene_symbols)

print("Count matrix successfully created")
print(dim(cts))

cell_types <- colData(patient_sce)[,"cell_type", drop=FALSE] |> 
    as.data.frame() 
print(dim(cell_types))

normal_cell_types <- (sce_normal$cell_type)[1]

print("Creating infercnv object")
infer <- CreateInfercnvObject(cts, 
    snakemake@input[["gene_order_file"]], 
    cell_types, 
    ref_group_names = normal_cell_types)

print("Running infercnv")
infercnv_obj <- infercnv::run(infer, 
    cutoff = snakemake@params[['cutoff']], 
    out_dir = snakemake@params[['infercnv_outdir']], 
    cluster_by_groups = TRUE, 
    analysis_mode = snakemake@params[['cnv_analysis_mode']], 
    denoise = snakemake@params[['denoise']], 
    sd_amplifier = ifelse(snakemake@params[['noise_logistic']], 3, 1.5), 
    noise_logistic = snakemake@params[['noise_logistic']], 
    tumor_subcluster_pval = 0.01,
    HMM = TRUE, 
    HMM_report_by = c('subcluster'), 
    HMM_type = 'i6', 
    BayesMaxPNormal = 0.5,
    no_prelim_plot = TRUE, 
    plot_steps = FALSE, 
    num_threads = as.numeric(snakemake@threads),
    leiden_resolution = snakemake@params[['leiden_res']], 
    png_res = 360)