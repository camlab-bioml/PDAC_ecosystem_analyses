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
  distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames("gene")

top.genes <- lapply(w, function(sig) {
  rownames(w)[sort(sig, index.return = TRUE, decreasing = TRUE)$ix[1:ntop]]
})
names(top.genes) <- names(w)

write_tsv(as.data.frame(top.genes), file = snakemake@output[["sig_top_gene_tsv"]])

top.genes <- Reduce(union, top.genes)

w.top.genes <- as.matrix(w)[top.genes,] %>%
  as.data.frame() %>%
  rownames_to_column("gene")

write_tsv(w.top.genes, file = snakemake@output[['sig_top_gene_loading_mtx']])

