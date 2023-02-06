suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(readr)
  library(sjstats)
  library(tidyr)
})

sigtopfreq.df <- read_tsv(snakemake@input[['sig_top_freq']])
sigloadings.df <- read_tsv(snakemake@input[['sig_loading_patient_summary']])
sigquantilepass.df <- read_tsv(snakemake@input[['sig_activation_freq']])

condition <- snakemake@wildcards[['condition']]
if (condition == "validated") condition = snakemake@wildcards[['subtype']]
if (condition == "collapsed") condition = paste(snakemake@wildcards[['subtype']], "Rep", sep = " ")

sigtopfreq.profiles <- sigtopfreq.df %>% 
  select(cohort, sample, contains(condition))

sigloadingmean.profiles <- sigloadings.df %>%
  select(-median) %>%
  pivot_wider(names_from = signature, values_from = mean) %>% 
  select(cohort, sample, contains(condition))

sigloadingmedian.profiles <- sigloadings.df %>%
  select(-mean) %>%
  pivot_wider(names_from = signature, values_from = median) %>% 
  select(cohort, sample, contains(condition))

sigquantilepass.profiles <- sigquantilepass.df %>%
  pivot_wider(names_from = signature, values_from = freq) %>%
  mutate_at(vars(contains(condition)), ~replace(., is.na(.), 0)) %>%
  select(cohort, sample, contains(condition)) %>%
  group_by(cohort, sample) %>%
  summarise_each(funs(sum)) %>%
  ungroup()

write_tsv(sigtopfreq.profiles, file = snakemake@output[['sig_top_freq_profiles']])
write_tsv(sigloadingmean.profiles, file = snakemake@output[['sig_loading_mean_profiles']])
write_tsv(sigloadingmedian.profiles, file = snakemake@output[['sig_loading_median_profiles']])
write_tsv(sigquantilepass.profiles, file = snakemake@output[['sig_act_freq_profiles']])





