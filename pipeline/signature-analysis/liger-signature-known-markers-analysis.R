suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(GeneOverlap)
})

w <- read_tsv(snakemake@input[['gene_loading_mtx']])
markerlists <- read_csv(snakemake@input[['marker_list']])

# process signature gene loading matrix
w$gene <- str_split(w$gene, pattern = "_", simplify = T)[,1]

print(paste0("number of unique genes in ", snakemake@wildcards[['condition']], " ", snakemake@wildcards[['subtype']], " signatures:"))
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


# Hypergeometric test ----------------------------------------------------------
## load signature gene loadings and known marker lists
celltype = snakemake@wildcards[['subtype']]
num_top_genes = snakemake@params[['num_top_genes_for_overlap_test']]

sig.gene.loading <- read_tsv(snakemake@input[['gene_loading_mtx']])
known.marker.lists <- read_csv(snakemake@input[['marker_list']])

## for collapsed version of signatures
if (snakemake@wildcards[['condition']] == "collapsed") celltype <- paste(celltype, "Rep", sep = " ")
if (snakemake@wildcards[['condition']] == "collapsed-scored-validation") celltype <- paste(celltype, "RepVal", sep = " ")

## make named lists
sig.top.genes <- lapply(seq(length(sig.gene.loading) - 1), function(sig) {
  signame <- paste(celltype, sig, sep = " ")
  genelist <- sig.gene.loading %>% 
    select(all_of(signame), gene) %>% 
    slice_max(order_by = get(signame), n = num_top_genes) %>% 
    pull(gene)
  str_split(genelist, "_", simplify = T)[,1]
})
names(sig.top.genes) <- names(sig.gene.loading)[seq(length(sig.gene.loading)-1)]

marker.lists <- lapply(known.marker.lists, function(markerlist) {
  markerlist[!is.na(markerlist)]
})

## construct GOM
gom.obj <- newGOM(gsetA = sig.top.genes, 
                  gsetB = marker.lists,
                  spec = "hg19.gene")

## save GOM
saveRDS(gom.obj, file = snakemake@output[['sig_gom']])










