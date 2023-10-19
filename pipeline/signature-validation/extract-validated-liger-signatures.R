suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(readr)
})

validation.metric = snakemake@params[['validation_metric']]
celltype <- snakemake@wildcards[['subtype']]

sigloading.df <- read_tsv(snakemake@input[['sig_loading_mtx']])
geneloading.df <- read_tsv(snakemake@input[['gene_loading_mtx']])
validated.sig.df <- read_tsv(snakemake@input[['validated_sig_df']])

validated.sig.df <- validated.sig.df %>%
  select(discovery, paste("validated", validation.metric, "sig", sep = ".")) %>%
  drop_na()
#validated.sig.df[is.na(validated.sig.df)] <- "invalid"

print(paste("Number of validated", celltype, "signatures:", nrow(validated.sig.df)))

# get validated signatures in signature loading matrix
names(sigloading.df) <- plyr::mapvalues(names(sigloading.df), 
                                        from = validated.sig.df[['discovery']], 
                                        to = validated.sig.df[[paste("validated", validation.metric, "sig", sep = ".")]])

sigloading.df <- sigloading.df %>%
  select(!contains("discovery"))

print(paste("Number of validated", celltype, "signatures in signature loading matrix:", ncol(sigloading.df %>% select(contains(celltype)))))

# get validated signatures in gene loading matrix
names(geneloading.df) <- plyr::mapvalues(names(geneloading.df),
                                         from = validated.sig.df[['discovery']],
                                         to = validated.sig.df[[paste("validated", validation.metric, "sig", sep = ".")]])

geneloading.df <- geneloading.df %>%
  select(!contains("discovery") & !contains("validation"))

print(paste("Number of validated", celltype, "signatures in gene loading matrix:", ncol(geneloading.df %>% select(contains(celltype)))))

# save the validated signatures
write_tsv(sigloading.df, file = snakemake@output[['validated_sig_loading_mtx']])
write_tsv(geneloading.df, file = snakemake@output[['validated_gene_loading_mtx']])


