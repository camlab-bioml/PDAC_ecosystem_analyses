suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(readr)
  library(BiocParallel)
  library(sjstats)
  library(dittoSeq)
  library(ggplotify)
  library(patchwork)
  library(cowplot)
  library(circlize)
})

celltypes <- snakemake@params[['celltypes']]
print(celltypes)

# load signature loading matrices
sig.loading.mtx.list <- lapply(celltypes, function(ct) {
  sig.loading.mtx.path <- grep(ct, snakemake@input[['sig_loading_mtx']], value = TRUE)[1]
  df <- read_tsv(sig.loading.mtx.path)
  names(df) <- gsub("Rep|RepVal", "Sig", names(df))
  names(df) <- gsub("discovery", paste0(ct, " discovery Sig"), names(df))
  names(df) <- gsub("validation", paste0(ct, " validation Sig"), names(df))
  df
})
names(sig.loading.mtx.list) <- celltypes

sig.loading.mtx.metadata.list <- lapply(sig.loading.mtx.list, function(df) {
  df |> select(where(is.character))
})

saveRDS(sig.loading.mtx.metadata.list, snakemake@output[['sig_loading_mtx_metadata']])

# summarize on sample level
## scale on cell level/restructure dataframes
sig.loading.mtx.list <- lapply(celltypes, function(ct) {
  sig.loading.mtx.df <- sig.loading.mtx.list[[ct]]
  sig.loading.mtx.metadata.df <- sig.loading.mtx.metadata.list[[ct]]
  cbind(sig.loading.mtx.df |> reframe(across(where(is.numeric), scale)),
        sig.loading.mtx.metadata.df)
})
names(sig.loading.mtx.list) <- celltypes

saveRDS(sig.loading.mtx.list, snakemake@output[['sig_loading_mtx_scaled']])

## get within sample variances
sig.loading.var.list <- lapply(celltypes, function(ct) {
  sig.loading.mtx.df <- sig.loading.mtx.list[[ct]]
  sig.loading.mtx.metadata.df <- sig.loading.mtx.metadata.list[[ct]]
  sig.loading.mtx.df |>
    reframe(across(where(is.numeric), var), .by = sample) |>
    mutate(cohort = plyr::mapvalues(sample,
                                    from = sig.loading.mtx.metadata.df$sample,
                                    to = sig.loading.mtx.metadata.df$cohort,
                        	    warn_missing = FALSE))
})
names(sig.loading.var.list) <- celltypes

sig.loading.var <- Reduce(cbind, lapply(sig.loading.var.list, function(df) {df |> summarise(across(where(is.numeric), median))})) |>
  as.matrix() |>
  t() |>
  set_colnames("Median of within sample var.") |>
  as.data.frame()

sig.loading.var.mtx <- Reduce(function(d1, d2) full_join(d1, d2, by = c("cohort", "sample")), sig.loading.var.list) |>
  as.data.frame()

print("sig.loading.var: ")
print(head(sig.loading.var))
print("sig.loading.var.mtx: ")
print(head(sig.loading.var.mtx))

saveRDS(sig.loading.var.mtx, snakemake@output[['sig_loading_var_mtx']])

## get sample means/medians
sig.loading.mean.list <- lapply(celltypes, function(ct) {
  sig.loading.mtx.df <- sig.loading.mtx.list[[ct]]
  sig.loading.mtx.metadata.df <- sig.loading.mtx.metadata.list[[ct]]
  sig.loading.mtx.df |> 
    reframe(across(where(is.numeric), mean), .by = sample) |>
    mutate(cohort = plyr::mapvalues(sample, 
                                    from = sig.loading.mtx.metadata.df$sample, 
                                    to = sig.loading.mtx.metadata.df$cohort,
                                    warn_missing = FALSE))
})
names(sig.loading.mean.list) <- celltypes

sig.loading.median.list <- lapply(celltypes, function(ct) {
  sig.loading.mtx.df <- sig.loading.mtx.list[[ct]]
  sig.loading.mtx.metadata.df <- sig.loading.mtx.metadata.list[[ct]]
  sig.loading.mtx.df |> 
    reframe(across(where(is.numeric), median), .by = sample) |>
    mutate(cohort = plyr::mapvalues(sample, 
                                    from = sig.loading.mtx.metadata.df$sample, 
                                    to = sig.loading.mtx.metadata.df$cohort,
                                    warn_missing = FALSE))
})
names(sig.loading.median.list) <- celltypes

## get variance of sample means/medians
sig.loading.mean.var.list <- lapply(celltypes, function(ct) {
  sig.loading.mean.df <- sig.loading.mean.list[[ct]]
  sig.loading.mtx.metadata.df <- sig.loading.mtx.metadata.list[[ct]]
  sig.loading.mean.df |> 
    reframe(across(where(is.numeric), var))
})
names(sig.loading.mean.var.list) <- celltypes

sig.loading.median.var.list <- lapply(celltypes, function(ct) {
  sig.loading.median.df <- sig.loading.median.list[[ct]]
  sig.loading.mtx.metadata.df <- sig.loading.mtx.metadata.list[[ct]]
  sig.loading.median.df |> 
    reframe(across(where(is.numeric), var))
})
names(sig.loading.median.var.list) <- celltypes


sig.loading.var <- sig.loading.var |>
  mutate(`Var. of sample means` = unlist(sig.loading.mean.var.list),
         `Var. of sample medians` = unlist(sig.loading.median.var.list)) |>
  rownames_to_column("signature")

print("sig.loading.var: ")
print(head(sig.loading.var))

write_tsv(sig.loading.var, snakemake@output[['sig_loading_var']])

## get between number of cells persample
sig.loading.ncell.list <- lapply(celltypes, function(ct) {
  sig.loading.mtx.df <- sig.loading.mtx.list[[ct]]
  sig.loading.mtx.metadata.df <- sig.loading.mtx.metadata.list[[ct]]
  sig.loading.mtx.df |> 
    select(!where(is.numeric), -cell_id) |>
    add_count(sample) |>
    distinct()
})
names(sig.loading.ncell.list) <- celltypes

saveRDS(sig.loading.ncell.list, snakemake@output[['sig_loading_ncell']])

## fill missing samples for some cell types with NA
sig.loading.median.mtx.df <- Reduce(function(d1, d2) full_join(d1, d2, by = c("cohort", "sample")), sig.loading.median.list)
sig.loading.mean.mtx.df <- Reduce(function(d1, d2) full_join(d1, d2, by = c("cohort", "sample")), sig.loading.mean.list)

write_tsv(sig.loading.median.mtx.df, snakemake@output[['sig_loading_median_mtx']])
write_tsv(sig.loading.mean.mtx.df, snakemake@output[['sig_loading_mean_mtx']])