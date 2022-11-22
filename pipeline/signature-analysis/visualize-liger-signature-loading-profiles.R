suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(readr)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
  library(corrplot)
  library(patchwork)
})

sig.profiles <- read_tsv(snakemake@input[['sig_profiles']])
condition <- snakemake@wildcards[['condition']]
if (condition == "validated") condition = snakemake@wildcards[['subtype']]
profile.flavour <- snakemake@wildcards[['profile']]

# construct the profile heatmap
col_ha <- columnAnnotation(Cohort = sig.profiles$cohort)

p <- Heatmap(sig.profiles %>% select(contains(condition)) %>% as.matrix() %>% t(), 
             column_labels = sig.profiles$sample, 
             name = profile.flavour,
             top_annotation = col_ha)

png(filename = snakemake@output[['sig_profiles_plot']], width = 10, height = 6, units = "in", res = 300)
draw(p, merge_legends = T)
dev.off()

png(filename = snakemake@output[['sig_profiles_corrplot']], width = 7, height = 7, units = "in", res = 300)
corrplot(cor(sig.profiles %>% select(contains(condition)) %>% select(order(names(.))) %>% as.matrix()))
dev.off()

