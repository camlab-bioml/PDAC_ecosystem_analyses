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
  library(scales)
})

celltypes <- snakemake@params[['celltypes']]
print(celltypes)

# load gene loading matrices
gene.loading.df.list <- lapply(celltypes, function(ct) {
  gene.loading.mtx.path <- grep(ct, snakemake@input[['gene_loading_mtx']], value = TRUE)[1]
  df <- read_tsv(gene.loading.mtx.path)
  df <- df |> mutate(gene = str_split(gene, "_ENSG", simplify = TRUE)[,1])
  df
})
names(gene.loading.df.list) <- celltypes

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

# how many sigs are left after validation and merging?
patient.collapsed.profiles.df <- read_tsv(snakemake@input[['patient_profiles_collapsed']])
patient.discovery.profiles.df <- read_tsv(snakemake@input[['patient_profiles_discovery']])
patient.validated.profiles.df <- read_tsv(snakemake@input[['patient_profiles_validated']])

sig.number.list <- list(
  Discovery = patient.discovery.profiles.df,
  Validated = patient.validated.profiles.df,
  Collapsed = patient.collapsed.profiles.df
)

sig.number.list <- lapply(sig.number.list, function(df) {
  data.frame (signature = names(df)[3:length(df)]) |>
    mutate(cell_type = gsub(" [0-9]$| Rep [0-9]$", "", signature)) |> 
    mutate(cell_type = gsub(" [0-9][0-9]$| Rep [0-9][0-9]$", "" , cell_type)) |>
    dplyr::count(cell_type)
})

sig.number.df <- left_join(sig.number.list$Discovery, sig.number.list$Validated, by = "cell_type", suffix = c("_Discovery", "_Validated"))
sig.number.df <- left_join(sig.number.df, sig.number.list$Collapsed, by = "cell_type")
names(sig.number.df) <- plyr::mapvalues(names(sig.number.df), from = c("n"), to = c("n_Merged"))

sig.number.df <- sig.number.df |>
  pivot_longer(starts_with("n_"), names_to = "Condition", values_to = "Number") |>
  mutate(Condition = factor(gsub("n_", "", Condition), levels = c("Discovery", "Validated", "Merged")))

saveRDS(sig.number.df, snakemake@output[['sig_number']])

# top markers for specific signatures
genes.unwanted <- c("MALAT1", "HEAT1", "XIST")

gene.loading.df.list <- lapply(gene.loading.df.list, function(gene.loading.df) {
  gene.loading.df %>%
    distinct(gene, .keep_all = TRUE) %>%
    filter(!(gene %in% genes.unwanted)) %>%
    column_to_rownames("gene")
})

top.gene.list <- lapply(gene.loading.df.list, function(gene.loading.df) {
  lapply(gene.loading.df, function(sig) {
    rownames(gene.loading.df)[sort(sig, index.return = TRUE, decreasing = TRUE)$ix[1:as.numeric(snakemake@params[['num_top_genes_to_show']])]]
  })
})

gene.loading.top.gene.df.list <- lapply(celltypes, function(ct) {
  gene.loading.df <- gene.loading.df.list[[ct]]
  top.genes <- Reduce(union, top.gene.list[[ct]])
  # add KRAS, BRCA1, BRCA2 for pancreatic epithelial cells
  if (ct == "pancreatic epithelial cell") {
    top.genes <- union(top.genes, c("KRAS"))
  }
  
  as.matrix(gene.loading.df)[top.genes,]
})
names(gene.loading.top.gene.df.list) <- celltypes


gene.loading.top.gene.df.to.plot.list <- list() 

# select signatures from epithelial cells
gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.list$`pancreatic epithelial cell` |> t() |> rescale() |> t()
top.genes.epi <- top.gene.list$`pancreatic epithelial cell`

gene.loading.top.gene.df.to.plot <- 
  gene.loading.top.gene.df.to.plot[Reduce(union,
                                          list(top.genes.epi$`pancreatic epithelial cell Rep 4`,
                                               top.genes.epi$`pancreatic epithelial cell Rep 7`,
                                               top.genes.epi$`pancreatic epithelial cell Rep 10`,
                                               top.genes.epi$`pancreatic epithelial cell Rep 11`,
                                               top.genes.epi$`pancreatic epithelial cell Rep 12`,
                                               #top.genes.epi$`pancreatic epithelial cell Rep 15`
                                               c("KRAS"))),
                                   c(4, 7, 10, 11, 12)] #, 15)]

gene.loading.top.gene.df.to.plot.list[['pancreatic epithelial cell']] <- gene.loading.top.gene.df.to.plot

# select signatures from fibroblasts
gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.list$fibroblast |> t() |> rescale() |> t()
top.genes.fibro <- top.gene.list$fibroblast

gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.to.plot[Reduce(union, 
                                                                            list(top.genes.fibro$`fibroblast Rep 2`,
                                                                                 top.genes.fibro$`fibroblast Rep 4`,
                                                                                 top.genes.fibro$`fibroblast Rep 6`,
                                                                                 top.genes.fibro$`fibroblast Rep 8`)),
                                                                     c(2, 4, 6, 8)]

gene.loading.top.gene.df.to.plot.list[['fibroblast']] <- gene.loading.top.gene.df.to.plot

# select signatures from CD8 T/NK cells
gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.list$`CD8-positive, alpha-beta T cell` |> t() |> rescale() |> t()
top.genes.cd8 <- top.gene.list$`CD8-positive, alpha-beta T cell`

gene.loading.top.gene.df.to.plot <- 
  gene.loading.top.gene.df.to.plot[Reduce(union, 
                                          list(top.genes.cd8$`CD8-positive, alpha-beta T cell Rep 2`,
                                               top.genes.cd8$`CD8-positive, alpha-beta T cell Rep 4`,
                                               top.genes.cd8$`CD8-positive, alpha-beta T cell Rep 7`,
                                               top.genes.cd8$`CD8-positive, alpha-beta T cell Rep 9`,
                                               top.genes.cd8$`CD8-positive, alpha-beta T cell Rep 10`)),
                                   c(2, 4, 7, 9, 10)]

gene.loading.top.gene.df.to.plot.list[['CD8-positive, alpha-beta T cell']] <- gene.loading.top.gene.df.to.plot

saveRDS(gene.loading.top.gene.df.to.plot.list, snakemake@output[['selected_sig_top_gene_loading_mtx_list']])