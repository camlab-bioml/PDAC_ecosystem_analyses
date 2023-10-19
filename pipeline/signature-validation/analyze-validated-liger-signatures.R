suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(readr)
})

celltype <- snakemake@wildcards[["subtype"]]

sigloading.df <- read_tsv(snakemake@input[['sig_loading_mtx']])
geneloading.df <- read_tsv(snakemake@input[['gene_loading_mtx']])

# compute validated signature signature loading correlations
sigloading.corr <- cor(as.matrix(sigloading.df |> select(contains(celltype))))

# compute validated signature gene loading correlations
geneloading.corr <- cor(as.matrix(geneloading.df |> select(contains(celltype))))
  
  
# save the validated signatures
write_tsv(as.data.frame(sigloading.corr), file = snakemake@output[['validated_sig_loading_corr']])
write_tsv(as.data.frame(geneloading.corr), file = snakemake@output[['validated_gene_loading_corr']])


