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

cohort.pal <- readRDS(snakemake@input[["cohort_pal"]])
cell_type_rename <- read_csv(snakemake@input[['cell_type_rename']])

sig.profiles <- read_tsv(snakemake@input[['sig_profiles']])
condition <- snakemake@wildcards[['condition']]
if (condition == "validated") condition = snakemake@wildcards[['subtype']]
if (condition == "collapsed") condition = paste(snakemake@wildcards[['subtype']], "Rep", sep = " ")
if (condition == "collapsed-scored-validation") condition = paste(snakemake@wildcards[['subtype']], "RepVal", sep = " ")

profile.flavour <- snakemake@wildcards[['profile']]

# construct the profile heatmap
col_ha <- columnAnnotation(Cohort = sig.profiles$cohort, col = list(Cohort = cohort.pal), show_annotation_name = FALSE)

mtx <- sig.profiles %>% select(contains(condition)) %>% as.matrix() %>% scale() %>% t()
mtx[is.na(mtx)] = 0

p <- Heatmap(mtx, 
             column_labels = sig.profiles$sample, 
             name = profile.flavour,
             top_annotation = col_ha,
             col = viridisLite::viridis(100, option = "B"))

png(filename = snakemake@output[['sig_profiles_plot']], width = 10, height = 6, units = "in", res = 360)
draw(p, merge_legends = T)
dev.off()

png(filename = snakemake@output[['sig_profiles_corrplot']], width = 7, height = 7, units = "in", res = 360)
corrplot(cor(sig.profiles %>% select(contains(condition)) %>% select(order(names(.))) %>% as.matrix() %>% scale()))
dev.off()

