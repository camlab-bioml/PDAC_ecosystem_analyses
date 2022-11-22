suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(readr)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(patchwork)
})

quantiles.cohort.df <- read_tsv(snakemake@input[['sig_loading_quantiles_cohort']])
quantiles.sample.df <- read_tsv(snakemake@input[['sig_loading_quantiles_sample']])
quantiles.signature.df <- read_tsv(snakemake@input[['sig_loading_quantiles_signature']])
quantiles.signature.cohort.df <- read_tsv(snakemake@input[['sig_loading_quantiles_signature_cohort']])
quantiles.signature.sample.df <- read_tsv(snakemake@input[['sig_loading_quantiles_signature_sample']])

# plot quantiles for each cohort
quantiles.cohort.df$probs <- as.character(quantiles.cohort.df$probs)

p9 <- ggplot(quantiles.cohort.df, aes(x = cohort, y = loading, fill = probs)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(y = "Loading", x = "Cohort", fill = "Quantile") + 
  theme_pubr()

# plot quantiles for each sample
quantiles.sample.df$probs <- as.character(quantiles.sample.df$probs)

p10 <- ggplot(quantiles.sample.df, aes(x = sample, y = loading, fill = probs)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(y = "Loading", x = "Sample", fill = "Quantile") + 
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# plot quantiles for each signature
quantiles.signature.df$probs <- as.character(quantiles.signature.df$probs)

p11 <- ggplot(quantiles.signature.df, aes(x = signature, y = loading, fill = probs)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(y = "Loading", x = "Signature", fill = "Quantile") + 
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 30, hjust=1))

# save the single variable bar plots
png(filename = snakemake@output[['sig_loading_quantiles_bar_plot']], width = 15, height = 8, units = "in", res = 300)
((p9 + p11) / p10) + 
  plot_layout(guides = "collect") & theme(legend.position = "right")
dev.off()

# plot quantiles for each signature and cohort
quantiles.signature.cohort.df$probs <- as.character(quantiles.signature.cohort.df$probs)

p12 <- ggplot(quantiles.signature.cohort.df, aes(x = signature, y = loading, color = probs)) + 
  geom_boxplot() + 
  labs(x = "Signature", y = "Signature loading", col = "Quantile") + 
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 30, hjust=1))

# plot quantiles for each signature and sample
quantiles.signature.sample.df$probs <- as.character(quantiles.signature.sample.df$probs)

p13 <- ggplot(quantiles.signature.sample.df, aes(x = signature, y = loading, color = probs)) + 
  geom_boxplot() + 
  labs(x = "Signature", y = "Signature loading", col = "Quantile") + 
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 30, hjust=1))

# save the multi-variable box plots
png(filename = snakemake@output[['sig_loading_quantiles_box_plot']], width = 15, height = 7, units = "in", res = 300)
(p12 + p13) + 
  plot_layout(guides = "collect") & theme(legend.position = "top")
dev.off()



