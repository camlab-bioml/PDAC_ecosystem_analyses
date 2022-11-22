suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(readr)
})

H.norm.df <- read_tsv(snakemake@input[['sig_loading_mtx']])

condition <- snakemake@wildcards[['condition']] #%>% stringr::str_to_title()
if (condition == "validated") condition = snakemake@wildcards[['subtype']]

sig_activation_threshold_quantile <- snakemake@params[['sig_activation_threshold_quantile']]

# get signature loading matrix
rownames(H.norm.df) <- H.norm.df$cell_id

H.norm.df <- H.norm.df %>% pivot_longer(cols = contains(condition), names_to = "signature", values_to = "loading")

write_tsv(H.norm.df, file = snakemake@output[['sig_loading_mtx_long']])

# analyze signature loading uniqueness
siguniqueness.df <- H.norm.df %>%
  group_by(cell_id) %>%
  mutate(loading_scaled = loading / max(loading)) %>%
  mutate(auc = sum(loading_scaled)) %>%
  mutate(auc_top2 = sum(sort(loading_scaled, decreasing = T)[1:2])) %>%
  mutate(top = signature[loading == max(loading)]) %>%
  mutate(second = ifelse(length(signature[loading == sort(loading, decreasing = T)[2]]) == 1, 
                         signature[loading == sort(loading, decreasing = T)[2]], 
                         NA))

write_tsv(siguniqueness.df, file = snakemake@output[['sig_loading_uniqueness']])

# check highest loaded signature in each cell
sigmax2.df <- H.norm.df %>%
  group_by(cell_id) %>%
  filter(loading == max(loading) | loading == sort(loading, decreasing = T)[2]) %>%
  mutate(rank = ifelse(loading == max(loading), "First", "Second")) %>%
  mutate(gap = max(loading) - min(loading)) %>%
  ungroup()

write_tsv(sigmax2.df, file = snakemake@output[['sig_loading_top_two']])

# check number of top signatures
sigtop.df <- sigmax2.df %>%
  dplyr::count(cohort, sample, signature, rank)

write_tsv(sigtop.df, file = snakemake@output[['sig_top_two_count']])

# check frequency of top signatures in each sample
sigtop.df <- sigmax2.df %>%
  filter(rank == "First") %>%
  group_by(cohort, sample, signature) %>%
  dplyr::summarise(n_cells = n()) %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  pivot_wider(names_from = signature, values_from = freq) %>%
  mutate_at(vars(contains(condition)), ~replace(., is.na(.), 0)) %>%
  summarise_each(funs(sum)) %>%
  ungroup()

write_tsv(sigtop.df, file = snakemake@output[['sig_top_freq']])

# summarize signature loadings in each patient
sigloadings.df <- H.norm.df %>%
  group_by(cohort, sample, signature) %>%
  summarise(mean = mean(loading), 
            median = median(loading)) %>%
  ungroup()

write_tsv(sigloadings.df, file = snakemake@output[['sig_loading_patient_summary']])

# get frequencies of activated signatures in each sample
quantile_probs <- c(0, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1)
my_quantile <- function(x, probs) {
  tibble(loading = quantile(x, probs), probs = probs)
}

quantiles.df <- H.norm.df %>% 
  group_by(signature) %>%
  summarise(my_quantile(loading, quantile_probs))
quantiles.df$probs <- as.character(quantiles.df$probs)
sigquantiles.df <- quantiles.df

sigquantilesthres.df <- sigquantiles.df %>% filter(probs == sig_activation_threshold_quantile)

sigquantilepass.df <- H.norm.df %>%
  group_by(cell_id) %>%
  mutate(passquantile = loading >= sigquantilesthres.df$loading[which(sigquantilesthres.df$signature == signature)]) %>%
  ungroup()

sigquantilepass.df <- sigquantilepass.df %>%
  select(-c(cell_id, loading)) %>%
  group_by(cohort, sample, signature) %>%
  dplyr::count(passquantile) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  filter(passquantile == T) %>%
  mutate(quantile = sig_activation_threshold_quantile)

write_tsv(sigquantilepass.df, file = snakemake@output[['sig_activation_freq']])


