suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(readr)
})

H.norm.df <- read_tsv(snakemake@input[['sig_loading_mtx_long']])

# check signature loading quantiles
quantile_probs <- c(0, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1)

my_quantile <- function(x, probs) {
  tibble(loading = quantile(x, probs), probs = probs)
}

## compute quantiles for each cohort
quantiles.df <- H.norm.df %>% 
  group_by(cohort) %>%
  summarise(my_quantile(loading, quantile_probs))

quantiles.df$probs <- as.character(quantiles.df$probs)

write_tsv(quantiles.df, file = snakemake@output[['sig_loading_quantiles_cohort']])

## compute quantiles for each sample
quantiles.df <- H.norm.df %>% 
  group_by(sample) %>%
  summarise(my_quantile(loading, quantile_probs))

quantiles.df$probs <- as.character(quantiles.df$probs)

write_tsv(quantiles.df, file = snakemake@output[['sig_loading_quantiles_sample']])

## compute quantiles for each signature
quantiles.df <- H.norm.df %>% 
  group_by(signature) %>%
  summarise(my_quantile(loading, quantile_probs))

quantiles.df$probs <- as.character(quantiles.df$probs)

sigquantiles.df <- quantiles.df

write_tsv(quantiles.df, file = snakemake@output[['sig_loading_quantiles_signature']])

## compute quantiles for each signature and cohort
quantiles.df <- H.norm.df %>% 
  group_by(signature, cohort) %>%
  summarise(my_quantile(loading, quantile_probs))

quantiles.df$probs <- as.character(quantiles.df$probs)

write_tsv(quantiles.df, file = snakemake@output[['sig_loading_quantiles_signature_cohort']])

## compute quantiles for each signature and sample
quantiles.df <- H.norm.df %>% 
  group_by(signature, sample) %>%
  summarise(my_quantile(loading, quantile_probs))

quantiles.df$probs <- as.character(quantiles.df$probs)

write_tsv(quantiles.df, file = snakemake@output[['sig_loading_quantiles_signature_sample']])



