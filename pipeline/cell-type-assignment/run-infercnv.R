suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(infercnv)
  library(Matrix)
  library(here)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(scales)
  library(bluster)
})

sce.tumor <- readRDS(snakemake@input[['sce_tumor']])
sce.normal <- readRDS(snakemake@input[['sce_normal']])
out.dir <- snakemake@params[['infercnv_outdir']]

# label cells
sce.normal$cell_type <- "Normal"
sce.tumor$cell_type <- colData(sce.tumor)[[paste0("predicted.", snakemake@params[['annot_level']])]]

# combine some cell type labels
sce.tumor$cell_type <- plyr::mapvalues(sce.tumor$cell_type, 
                                       from = c("alpha", "beta", "delta", "gamma", "activated_stellate", "quiescent_stellate"),
                                       to = c("endocrine", "endocrine", "endocrine", "endocrine", "stellate", "stellate"))

# tidy up metadata for cells
common.coldata <- intersect(names(colData(sce.normal)), names(colData(sce.tumor)))
colData(sce.normal) <- colData(sce.normal)[,common.coldata]
colData(sce.tumor) <- colData(sce.tumor)[,common.coldata]

reducedDims(sce.normal) <- NULL
reducedDims(sce.tumor) <- NULL

# tidy up genes for reference and observation 
rownames(sce.normal) <- paste(rowData(sce.normal)[['Symbol']], rowData(sce.normal)[['ensembl_id']], sep = "_")
rownames(sce.tumor) <- paste(rowData(sce.tumor)[['Symbol']], rowData(sce.tumor)[['ensembl_id']], sep = "_")

common.gene <- intersect(rownames(sce.normal), rownames(sce.tumor))
sce.normal <- sce.normal[common.gene,]
sce.tumor <- sce.tumor[common.gene,]

rownames(sce.normal) <- str_split(rownames(sce.normal), pattern = "_", simplify = T)[,1]
rownames(sce.tumor) <- str_split(rownames(sce.tumor), pattern = "_", simplify = T)[,1]

# combine SCEs
sce.combined <- cbind(sce.tumor, sce.normal)

# get the raw count matrix
counts.mtx <- counts(sce.combined)

# save the tumor annotation file
write.table(colData(sce.combined)[,'cell_type', drop=F], file = paste0(out.dir, 'tumor_classes.txt'), sep="\t", col.names=FALSE, quote = FALSE)

# save the gene order file
print(snakemake@input[['gene_order_file']])
order_csv <-  read.csv(snakemake@input[['gene_order_file']], sep = '\t', header=FALSE)
order_csv$V1 <- sapply(lapply(lapply(order_csv$V1, strsplit, "|", fixed=TRUE), '[[',1),'[',1)
order_csv <- order_csv[!duplicated(order_csv$V1),,drop=FALSE]

write.table(order_csv, file = paste0(out.dir, 'gene_order.txt'), sep="\t", col.names=FALSE, row.names=FALSE, quote = FALSE)

# create inferCNV object
infer <- CreateInfercnvObject(
  raw_counts_matrix= counts.mtx,
  delim="\t",
  annotations_file= paste0(out.dir, 'tumor_classes.txt'),
  gene_order_file= paste0(out.dir, 'gene_order.txt'),
  ref_group_names= c("Normal"))

# run inferCNV
infercnv_obj <- infercnv::run(
  infercnv_obj = infer,
  cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  min_cells_per_gene = 3,
  window_length = 101,
  out_dir = out.dir,
  cluster_by_groups = T, 
  plot_steps = F,
  denoise = T,
  HMM = T,
  HMM_type = "i6",
  analysis_mode = "samples",
  tumor_subcluster_partition_method = "leiden",
  no_prelim_plot = T,
  num_threads = as.numeric(snakemake@threads), 
  png_res = 321,
  resume_mode = T,
  save_rds = F,
  save_final_rds = T
)


