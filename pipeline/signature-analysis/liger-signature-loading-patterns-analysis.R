suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(readr)
  library(sjstats)
  library(tidyr)
  library(singleCellTK)
})

# load signature loadings in single cells
sig.tops <- read_tsv(snakemake@input[["sig_loading_top_two"]])

celltype <- snakemake@wildcards[["subtype"]]
condition <- snakemake@wildcards[["condition"]]

# load sces
if (condition == "collapsed") sce <- readRDS(snakemake@input[["sce_dis"]])
if (condition == "collapsed-scored-validation") sce <- readRDS(snakemake@input[["sce_val"]])

dim.red.plot <- snakemake@params[["dim_red_plot"]]

print("Available reduced dimensionality reduction matrices for plotting: ")
reducedDimNames(sce)

#profile.flavour <- snakemake@wildcards[["profile"]]

# construct the pattern umap
sig.top.one <- sig.tops |> filter(rank == "First")
sig.top.two <- sig.tops |> filter(rank == "Second")

## get Seurat UMAP
print(paste0("Running Seurat UMAP for ", celltype, " under condition: ", condition))
sce <- runSeuratUMAP(sce, nNeighbors = 30, useReduction = "pca", seed = 1, verbose = TRUE)
print("Seurat UMAP completed")

# get desired metadata
df.redim <- data.frame(
  Cell_ID = colnames(sce),
  Sample = sce$sample,
  Cohort = sce$cohort,
  Cell_type = celltype,
  top_sig = plyr::mapvalues(colnames(sce),
    from = sig.top.one$cell_id,
    to = sig.top.one$signature,
    warn_missing = FALSE
  ),
  second_sig = plyr::mapvalues(colnames(sce),
    from = sig.top.two$cell_id,
    to = sig.top.two$signature,
    warn_missing = FALSE
  ),
  top_sig_loading = plyr::mapvalues(colnames(sce),
    from = sig.top.one$cell_id,
    to = sig.top.one$loading,
    warn_missing = FALSE
  ),
  second_sig_loading = plyr::mapvalues(colnames(sce),
    from = sig.top.two$cell_id,
    to = sig.top.two$loading,
    warn_missing = FALSE
  ),
  UMAP_1 = reducedDim(sce, dim.red.plot)[, 1],
  UMAP_2 = reducedDim(sce, dim.red.plot)[, 2]
)

# tidyup metadata
df.redim <- df.redim |>
  mutate(
    top_sig = ifelse(grepl(celltype, top_sig), gsub(" Rep | RepVal ", " ", top_sig), NA),
    second_sig = ifelse(grepl(celltype, second_sig), gsub(" Rep | RepVal ", " ", second_sig), NA),
    top_sig_loading = ifelse(grepl("_|-", top_sig_loading), NA, top_sig_loading),
    second_sig_loading = ifelse(grepl("_|-", second_sig_loading), NA, second_sig_loading)
  )

df.redim$top_sig_loading <- df.redim$top_sig_loading |> as.numeric()
df.redim$second_sig_loading <- df.redim$second_sig_loading |> as.numeric()

# save df.redim
write_tsv(df.redim, file = snakemake@output[["dimred_with_top_two_sig_loadings"]])
