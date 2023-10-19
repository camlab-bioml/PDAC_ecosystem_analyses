suppressPackageStartupMessages({
  library(SingleCellExperiment) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(magrittr)
  library(stats)
  library(ggpubr)
})

celltypes <- snakemake@params[["celltypes"]]
condition <- snakemake@wildcards[["condition"]]

print("Cell types used for Stan model:")
print(celltypes)
print("Condition used for Stan model:")
print(condition)

number.of.niches = snakemake@params[["number_of_niches"]] |> as.numeric()
min.number.of.cells.per.sample = snakemake@params[["min_number_of_cells_per_sample"]]
max.number.of.cells.per.sample = snakemake@params[["max_number_of_cells_per_sample"]]

print("Number of niches used in Stan model:")
print(number.of.niches)
print("Minimum number of cells per sample used in Stan model:")
print(min.number.of.cells.per.sample)
print("Maximum number of cells per sample used in Stan model:")
print(max.number.of.cells.per.sample)

# load data
df.sig.list <- lapply(snakemake@input[['celltype_sig_loading_mtxs']], read_tsv)
names(df.sig.list) <- celltypes

df.sig.mean <- read_tsv(snakemake@input[['sig_loading_profiles']])

# removing samples with too few cells
lapply(celltypes, function(ct) {
  holder <- df.sig.list[[ct]]
  p <- ggplot(holder, aes(x = sample, fill = cohort)) +
    geom_bar() + 
    scale_y_log10() + 
    labs(title = ct) + 
    theme_pubr(x.text.angle = 45)
  ggsave(filename = paste0("number-of-cells-per-sample-", ct, "-", condition, ".png"), path = snakemake@params[["prep_fig_path"]],
         plot = p, device = "png", bg = "white", width = 10, height = 5, units = "in", dpi = "retina")
})

df.sig.list <- lapply(df.sig.list, function(df.sig) {
  df.sig %>%
    group_by(sample) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    dplyr::filter(count >= min.number.of.cells.per.sample) %>%
    dplyr::select(!count)
})

samples.common <- Reduce(intersect, lapply(df.sig.list, function(df.sig) unique(df.sig$sample)))
samples.union <- Reduce(union, lapply(df.sig.list, function(df.sig) unique(df.sig$sample)))
samples.union

# construct elements of Stan data list
stan.data.list <- lapply(names(df.sig.list), function(ct) {
  df.sig <- df.sig.list[[ct]] %>% 
    #dplyr::filter(sample %in% samples.common) %>% 
    slice_sample(n = max.number.of.cells.per.sample, by = sample, replace = FALSE) %>%
    mutate(sample = factor(sample, levels = unique(sample))) %>% 
    arrange(sample)
  list(
    N = nrow(df.sig),
    P = length(unique(df.sig$sample)),
    K = length(grep(ct, names(df.sig), value = TRUE)),
    Y = df.sig %>% dplyr::select(contains(ct)) %>% as.matrix(),
    y = (df.sig %>% dplyr::select(contains(ct)) %>% as.matrix()) / sd(df.sig %>% dplyr::select(contains(ct)) %>% as.matrix()),
    x = plyr::mapvalues(df.sig$sample, 
                        from = samples.union,
                        to = seq_along(samples.union) %>% as.numeric(),
                        warn_missing = FALSE)
  )
})
names(stan.data.list) <- names(df.sig.list)

# build the stan data list for model fitting
stan.data.base <- list(
  C = length(stan.data.list),
  L = number.of.niches,
  P = length(samples.union)
)

stan.data <- stan.data.base

j.tracker = 1

for (ct in celltypes) {
  stan.data.ct <- list(
    j = j.tracker,
    N = stan.data.list[[ct]]$N,
    P = stan.data.list[[ct]]$P,
    K = stan.data.list[[ct]]$K,
    y = stan.data.list[[ct]]$y,
    x = stan.data.list[[ct]]$x
  )
  names(stan.data.ct) <- paste0(names(stan.data.ct), "_", ct)
  j.tracker = j.tracker + stan.data.ct$K
  
  stan.data <- c(stan.data, stan.data.ct)
}

rm(j.tracker, ct, stan.data.ct)

saveRDS(stan.data, file = snakemake@output[["stan_data"]])

# tidy up observed patient level signature loading mean mtrix
sigs.order <- lapply(stan.data.list, function(df) {colnames(df$y)}) %>% unlist()
sigs.order

df.sig.mean <- df.sig.mean %>%
  slice(order(factor(sample, levels = samples.union))) %>%
  column_to_rownames("sample") %>%
  select(all_of(sigs.order)) %>%
  mutate_if(is.numeric, scales::rescale) %>%
  replace(is.na(.), 0)

saveRDS(df.sig.mean, file = snakemake@output[["df_sig_mean"]])