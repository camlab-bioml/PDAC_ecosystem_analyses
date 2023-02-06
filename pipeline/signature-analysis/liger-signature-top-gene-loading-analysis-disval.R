suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
})

w <- read_tsv(snakemake@input[['gene_loading_mtx']])
ntop <- snakemake@params[['num_top_genes']] %>% as.numeric()

w$gene <- str_split(w$gene, pattern = "_", simplify = T)[,1]

print(paste0("number of unique genes in ", snakemake@wildcards[['subtype']], " signatures:"))
print(length(unique(w$gene)))

w <- w %>%
  distinct(gene, .keep_all = T)

gene.holder <- w$gene
w$gene <- NULL

rownames(w) <- gene.holder
rm(gene.holder)

top.genes <- lapply(w, function(sig) {
  rownames(w)[sort(sig, index.return=TRUE, decreasing=TRUE)$ix[1:ntop]]
})
top.genes <- Reduce(union, top.genes)

w.top.genes <- as.matrix(w)[top.genes,]
top.genes <- rownames(w.top.genes)

w.top.genes <- w.top.genes %>% as.data.frame()
w.top.genes$gene <- top.genes

write_tsv(w.top.genes, file = snakemake@output[['sig_top_gene_loading_mtx']])


