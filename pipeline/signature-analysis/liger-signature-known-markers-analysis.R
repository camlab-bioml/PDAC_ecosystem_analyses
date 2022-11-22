suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
})

w <- read_tsv(snakemake@input[['gene_loading_mtx']])
markerlists <- read_csv(snakemake@input[['marker_list']])

# process signature gene loading matrix
w$gene <- str_split(w$gene, pattern = "_", simplify = T)[,1]

print(paste0("number of unique genes in ", snakemake@wildcards[['subtype']], " signatures:"))
print(length(unique(w$gene)))

w <- w %>%
  distinct(gene, .keep_all = T)

gene.holder <- w$gene
w$gene <- NULL
rownames(w) <- gene.holder
rm(gene.holder)

# process marker lists to get marker group labels
row_split <- sapply(seq(length(markerlists)), function(i) {
  name <- names(markerlists)[i]
  markers <- markerlists[[i]][!is.na(markerlists[[i]])]
  rep(name, length(markers))
}) %>% Reduce(c, .)

# subset gene loading matrix to marker genes
w.markers <- lapply(markerlists, function(markers) {
  markers <- markers[!is.na(markers)]
  w.markers.mtx <- w[match(markers, rownames(w)),] %>% as.matrix()
  rownames(w.markers.mtx) <- markers
  w.markers.mtx
})
w.markers <- Reduce(rbind, w.markers)

w.markers <- w.markers %>% as.data.frame()
w.markers$split <- row_split
w.markers$gene <- rownames(w.markers)

# save the marker loading matrix
write_tsv(w.markers, file = snakemake@output[['sig_known_gene_loading_mtx']])

# =============================================================================

w <- read_tsv(snakemake@input[['gene_loading_mtx_validated']])
markerlists <- read_csv(snakemake@input[['marker_list']])

# process signature gene loading matrix
w$gene <- str_split(w$gene, pattern = "_", simplify = T)[,1]

print(paste0("number of unique genes in validated ", snakemake@wildcards[['subtype']], " signatures:"))
print(length(unique(w$gene)))

w <- w %>%
  distinct(gene, .keep_all = T)

gene.holder <- w$gene
w$gene <- NULL
rownames(w) <- gene.holder
rm(gene.holder)

# process marker lists to get marker group labels
row_split <- sapply(seq(length(markerlists)), function(i) {
  name <- names(markerlists)[i]
  markers <- markerlists[[i]][!is.na(markerlists[[i]])]
  rep(name, length(markers))
}) %>% Reduce(c, .)

# subset gene loading matrix to marker genes
w.markers <- lapply(markerlists, function(markers) {
  markers <- markers[!is.na(markers)]
  w.markers.mtx <- w[match(markers, rownames(w)),] %>% as.matrix()
  rownames(w.markers.mtx) <- markers
  w.markers.mtx
})
w.markers <- Reduce(rbind, w.markers)

w.markers <- w.markers %>% as.data.frame()
w.markers$split <- row_split
w.markers$gene <- rownames(w.markers)

# save the marker loading matrix
write_tsv(w.markers, file = snakemake@output[['sig_known_gene_loading_mtx_validated']])

