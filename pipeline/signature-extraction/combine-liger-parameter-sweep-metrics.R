suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(tidyr)
  library(stringi)
  library(readr)
})

tsvpathlist <- snakemake@input[['tsv']]

tsvlist <- lapply(tsvpathlist, read_tsv)

tsv_combined <- Reduce(bind_rows, tsvlist) %>%
  drop_na() %>%
  mutate(parameter = stri_split_fixed(as.character(parameter), "_")) %>%
  unnest(parameter) %>%
  mutate(parameter = stri_trim_both(parameter)) %>%
  separate(parameter, into = c("param", "val"), sep = ":")

params <- unique(tsv_combined$param)

tsv_combined <- tsv_combined %>%
  pivot_wider(names_from = param, values_from = val) %>%
  unnest(params)

write_tsv(tsv_combined, file = snakemake@output[['tsv']])