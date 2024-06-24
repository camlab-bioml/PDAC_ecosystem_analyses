suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(scuttle)
  library(readr)
})

w <- read_tsv(snakemake@input[['gene_loading_mtx']])
ntop <- snakemake@params[['num_top_genes']] %>% as.numeric()

w$gene <- str_split(w$gene, pattern = "_", simplify = T)[,1]

print(paste0("number of unique genes in ", snakemake@wildcards[['condition']], " ", snakemake@wildcards[['subtype']], " signatures:"))
print(length(unique(w$gene)))

w <- w %>%
  distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames("gene")

top.genes <- lapply(w, function(sig) {
  rownames(w)[sort(sig, index.return = TRUE, decreasing = TRUE)$ix[1:ntop]]
})
names(top.genes) <- names(w)

write_tsv(as.data.frame(top.genes), file = snakemake@output[['sig_top_gene_tsv']])

top.genes <- Reduce(union, top.genes)

w.top.genes <- as.matrix(w)[top.genes,] %>%
  as.data.frame() %>%
  rownames_to_column("gene")

write_tsv(w.top.genes, file = snakemake@output[['sig_top_gene_loading_mtx']])

# further analyze top genes loadings, ranks, within cell type expression, etc.
## load expression data
if (snakemake@wildcards[['condition']] == "collapsed-scored-validation") {
  sce <- readRDS(snakemake@input[['sce_val']])
} else {
  sce <- readRDS(snakemake@input[['sce_dis']])
}

## load signature loadings
sig.loading.df <- read_tsv(snakemake@input[['sig_loading_mtx']])

## analysis
w.top.genes.df <- w |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  # filter(gene %in% top.genes) |>
  pivot_longer(cols = -gene, names_to = "signature", values_to = "loading") |>
  mutate(celltype = str_split(signature, pattern = " Rep | RepVal | [0-9]$| [0-9][0-9]$", simplify = TRUE)[,1]) |>
  #mutate(celltype = str_replace_all(celltype, pattern = " ", replacement = "_")) |>
  # mutate(cell_type = plyr::mapvalues(cell_type,
  #                                    from = cell_type_rename$old_name,
  #                                    to = cell_type_rename$new_name)) |>
  group_by(signature) |>
  mutate(rank_loading = rank(-loading)) |>
  slice_max(order_by = loading, n = ntop) |>
  ungroup()
  # mutate(loading = rescale_1_99(loading)) |>
  # mutate(loading = ifelse(loading < 0, 0, loading)) |>
  # mutate(loading = ifelse(loading > 1, 1, loading)) |>
  # mutate(loading = round(loading, 2)) |>
  # mutate(loading = as.character(loading)) |>
  # mutate(loading = paste0(loading, "%"))

# add within cell type expression
length(unique(w.top.genes.df$gene))
length(intersect(unique(w.top.genes.df$gene), rownames(sce)))

genes.common <- intersect(unique(w.top.genes.df$gene), rownames(sce))
logcounts.mtx <- assay(sce, "logcounts")#[genes.common,]

within.celltype.exprs.df <- data.frame(
  gene = rownames(logcounts.mtx),
  within_celltype_exprs_mean = rowMeans(logcounts.mtx),
  within_celltype_exprs_mean_rank = rank(-rowMeans(logcounts.mtx)),
  within_celltype_exprs_var = rowVars(logcounts.mtx)
)

w.top.genes.df <- w.top.genes.df %>%
  left_join(within.celltype.exprs.df, by = "gene")

## add signature loadings correlation with within cell type expression
sig.loading.df <- sig.loading.df |> arrange(factor(cell_id, levels = colnames(logcounts.mtx)))

sig.loading.gene.exprs.cor.mtx <- cor(x = sig.loading.df |> select(contains(snakemake@wildcards[['subtype']])) |> as.matrix(),
                                      y = assay(sce, "logcounts") |> as.matrix() |> t(), 
                                      use = "pairwise.complete.obs")
dim(sig.loading.gene.exprs.cor.mtx)
sig.loading.gene.exprs.cor.mtx <- t(sig.loading.gene.exprs.cor.mtx)

sig.loading.gene.exprs.cor.df <- sig.loading.gene.exprs.cor.mtx |> as.data.frame() |> rownames_to_column("gene") |> 
  pivot_longer(cols = -gene, names_to = "signature", values_to = "within_celltype_exprs_sig_loading_corr") |> 
  # mutate(celltype = str_split(signature, pattern = " Rep | RepVal | [0-9]$| [0-9][0-9]$", simplify = TRUE)[,1]) |>
  group_by(signature) |>
  mutate(rank_exprs_sig_loading_corr = rank(-within_celltype_exprs_sig_loading_corr)) |>
  #slice_max(order_by = cor, n = ntop) |>
  ungroup()

w.top.genes.df <- w.top.genes.df |> left_join(sig.loading.gene.exprs.cor.df, by = c("gene", "signature"))

# rearrange order of columns
w.top.genes.df <- w.top.genes.df |> select(celltype, signature, gene, 
                                           loading, rank_loading, 
                                           within_celltype_exprs_sig_loading_corr, rank_exprs_sig_loading_corr, 
                                           within_celltype_exprs_mean, within_celltype_exprs_mean_rank, within_celltype_exprs_var)

write_tsv(w.top.genes.df, file = snakemake@output[['sig_top_gene_analysis_df']])