suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(readr)
  library(sjstats)
  library(tidyr)
  library(SingleCellExperiment)
  library(singleCellTK)
})

# load signature loadings in single cells
sig.tops <- read_tsv(snakemake@input[["sig_loading_top_two"]])

celltype <- snakemake@wildcards[["subtype"]]
condition <- snakemake@wildcards[["condition"]]

# load sces
if (condition == "collapsed") sce <- readRDS(snakemake@input[["sce_dis"]])
if (condition == "collapsed-scored-validation") sce <- readRDS(snakemake@input[["sce_val"]])

## remove zero expression genes
sce <- sce[Matrix::rowSums(assay(sce, "logcounts") > 0) > 0, ]


# get top signature in single cells
sig.top.one <- sig.tops |> filter(rank == "First")
#sig.top.two <- sig.tops |> filter(rank == "Second")

# get desired expression
exprs.mtx <- assay(sce, "logcounts") |> as.matrix() |> t()

# tidyup expression
exprs.mtx <- exprs.mtx |> as.data.frame() 

which(duplicated(names(exprs.mtx)))
names(exprs.mtx)[duplicated(names(exprs.mtx))]
# for B cells, remove the duplicated gene names
exprs.mtx <- exprs.mtx |> select(-names(exprs.mtx)[duplicated(names(exprs.mtx))])

exprs.mtx <- exprs.mtx  |>
  mutate(top_sig = plyr::mapvalues(rownames(exprs.mtx),
                                   from = sig.top.one$cell_id,
                                   to = sig.top.one$signature,
                                   warn_missing = FALSE)
        ) |>
  mutate(cohort = sce$cohort,
         sample = sce$sample) |>
  filter(grepl("Rep|RepVal", top_sig))

# average expression for all genes per top signature per sample per cohort
exprs.mtx <- exprs.mtx |> 
  group_by(top_sig) |>
  summarise_all(mean) |>
  ungroup()

# save top signture gene expressions in samples
write_tsv(exprs.mtx, file = snakemake@output[["sig_gene_expr_mtx"]])
