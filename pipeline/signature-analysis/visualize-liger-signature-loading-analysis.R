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

H.norm.df <- read_tsv(snakemake@input[['sig_loading_mtx_long']])
siguniqueness.df <- read_tsv(snakemake@input[['sig_loading_uniqueness']])
sigmax2.df <- read_tsv(snakemake@input[['sig_loading_top_two']])
sigtop.df <- read_tsv(snakemake@input[['sig_top_two_count']])
sigtopfreq.df <- read_tsv(snakemake@input[['sig_top_freq']])

condition <- snakemake@wildcards[['condition']] 
if (condition == "validated") condition = snakemake@wildcards[['subtype']]
if (condition == "collapsed") condition = paste(snakemake@wildcards[['subtype']], "Rep", sep = " ")

# make quantile-quantile plot for signature loading
p1 <- ggplot(H.norm.df, aes(sample = loading)) +
  stat_qq(pch = 1) +
  stat_qq_line(col = "steelblue", lwd = 1) +
  labs(x = "Theoretical Quantiles", y = "Signature loadings") + 
  theme_pubr()

p2 <- ggplot(H.norm.df, aes(sample = loading)) +
  stat_qq(pch = 1, aes(colour = signature)) +
  stat_qq_line(col = "grey30", lwd = 1) +
  labs(x = "Theoretical Quantiles", y = "Signature loadings",
       col = "Signature") + 
  theme_pubr(legend = "right")

png(filename = snakemake@output[['sig_loading_qqnorm_plot']], width = 10, height = 5, units = "in", res = 300)
p1 + p2 + 
  plot_annotation(title = paste0("Normal Q-Q plot for ", snakemake@wildcards[['subtype']], " signature loadings"))
dev.off()

# plot uniqueness of signatures
tmp.siguniqueness.df <- siguniqueness.df %>%
  select(-signature, -loading, -loading_scaled) %>%
  distinct()

p3 <- ggplot(tmp.siguniqueness.df, aes(x = auc, color = top, fill = top)) + 
  geom_density(alpha = .1) + 
  #geom_density(data = tmp.siguniqueness.df, mapping = aes(x = auc), alpha = .1, colour = "grey70") + 
  facet_grid(second ~ top) +
  labs(title = paste0("AUC grouped by top signature (column) and secondary signature (row) in ", snakemake@wildcards[['subtype']], " cells"),
       col = "Top signature", fill = "Top signature",
       x = "AUC", y = "Density") +
  theme_pubr()

p4 <- ggplot(tmp.siguniqueness.df, aes(x = auc, color = second, fill = second)) + 
  geom_density(alpha = .1) + 
  #geom_density(data = tmp.siguniqueness.df, mapping = aes(x = auc), alpha = .1, colour = "grey70") + 
  facet_wrap(~ top, scales = "free", ncol = 2) +
  labs(col = "second top signature", fill = "second top signature",
       x = "AUC", y = "Density",
       title = paste0("AUC grouped by top signature in the cell")) +
  theme_pubr()

p5 <- ggplot(tmp.siguniqueness.df, aes(x = top, y = auc, color = second)) + 
  geom_boxplot() + 
  #facet_wrap(~ top, scales = "free") +
  labs(col = "second top signature", 
       x = "Top signature",
       y = "AUC - all signatures") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

p6 <- ggplot(tmp.siguniqueness.df, aes(x = top, y = auc)) + 
  geom_boxplot() +
  labs(x = "Top signature",
       y = "AUC - all signatures") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

p7 <- ggplot(tmp.siguniqueness.df, aes(x = top, y = auc_top2, color = second)) + 
  geom_boxplot() + 
  #facet_wrap(~ top, scales = "free") +
  labs(col = "second top signature", 
       x = "Top signature",
       y = "AUC - top 2 signatures") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

p8 <- ggplot(tmp.siguniqueness.df, aes(x = top, y = auc_top2)) + 
  geom_boxplot() +
  labs(x = "Top signature",
       y = "AUC - top 2 signatures") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

png(filename = snakemake@output[['sig_loading_auc_density_grid_plot']], width = 10, height = 10, units = "in", res = 300)
p3
dev.off()

png(filename = snakemake@output[['sig_loading_auc_density_grouped_plot']], width = 10, height = 8, units = "in", res = 300)
p4
dev.off()

png(filename = snakemake@output[['sig_loading_auc_box_grouped_plot']], width = 10, height = 8, units = "in", res = 300)
(p5 / p7) + 
  plot_annotation(title = paste0("AUC for scaled ", snakemake@wildcards[['subtype']], " signature loading")) + 
  plot_layout(guides = "collect") & theme(legend.position = "right")
dev.off()

png(filename = snakemake@output[['sig_loading_auc_box_summary_plot']], width = 8, height = 8, units = "in", res = 300)
(p6 / p8) +
  plot_annotation(title = paste0("AUC for scaled ", snakemake@wildcards[['subtype']], " signature loading"))
dev.off()

rm(tmp.siguniqueness.df)

# plot loadings of highest and second highest loaded signatures in cells
p9 <- ggplot(sigmax2.df, aes(x = signature, y = loading, color = rank, fill = rank)) + 
  geom_boxplot(alpha = .3) + 
  #facet_wrap(~ cohort, nrow = 2) +
  labs(col = "Signature rank", fill = "Signature rank",
       y = "Loading", x = NULL) + 
  theme_pubr() + 
  theme(axis.text.x = element_blank())

p10 <- ggplot(sigmax2.df, aes(x = signature, y = gap, color = rank, fill = rank)) + 
  geom_boxplot(alpha = .3) + 
  #facet_wrap(~ cohort, nrow = 2) +
  labs(col = "Signature rank", fill = "Signature rank",
       y = "Loading gap", x = "Signature") + 
  theme_pubr(x.text.angle = 45)

p11 <- p9 + facet_wrap(~ cohort, nrow = 2)
p12 <- p10 + facet_wrap(~ cohort, nrow = 2)

png(filename = snakemake@output[['sig_top_loading_box_summary_plot']], width = 10, height = 8, units = "in", res = 300)
(p9 / p10) + 
  plot_annotation(title = paste0("Loading of top 2 ", snakemake@wildcards[['subtype']], " signatures")) + 
  plot_layout(guides = "collect") & theme(legend.position = "right")
dev.off()

png(filename = snakemake@output[['sig_top_loading_box_grouped_plot']], width = 10, height = 10, units = "in", res = 300)
(p11 / p12) + 
  plot_annotation(title = paste0("Loading of top 2 ", snakemake@wildcards[['subtype']], " signatures in each cohort")) + 
  plot_layout(guides = "collect") & theme(legend.position = "right")
dev.off()

# plot numbers of highest and second highest loaded signatures in samples
p13 <- ggplot(sigtop.df, aes(x = signature, y = n, color = rank, fill = rank)) + 
  geom_boxplot(alpha = .3) + 
  facet_wrap(~ cohort, nrow = 2, scales = "free") +
  labs(col = "Signature rank", fill = "Signature rank",
       y = "Number of cells", x = "Signature") + 
  theme_pubr(x.text.angle = 45)

p14 <- ggplot(sigtop.df, aes(x = rank, y = n, color = rank, fill = rank)) + 
  geom_boxplot(alpha = .3) + 
  facet_wrap(~ signature, scales = "free") +
  labs(col = "Signature rank", fill = "Signature rank",
       y = "Number of cells", x = "Signature") + 
  theme_pubr(x.text.angle = 0)

png(filename = snakemake@output[['sig_top_count_box_grouped_plot']], width = 10, height = 8, units = "in", res = 300)
p13 + 
  plot_annotation(title = paste0("Number of top 2 ", snakemake@wildcards[['subtype']], " signatures in each sample"))
dev.off()

png(filename = snakemake@output[['sig_top_count_box_summary_plot']], width = 10, height = 8, units = "in", res = 300)
p14 + 
  plot_annotation(title = paste0("Number of top 2 ", snakemake@wildcards[['subtype']], " signatures"))
dev.off()

# plot frequency of highest loaded signatures in samples
tmp.sigtopfreq.df <- sigtopfreq.df %>%
  pivot_longer(cols = contains(condition), names_to = "signature", values_to = "freq")

p15 <- ggplot(tmp.sigtopfreq.df, aes(x = sample, y = freq, fill = signature)) + 
  geom_bar(position = "stack", stat = "identity") + 
  #facet_wrap(~ cohort, scales = "free") +
  labs(fill = "Signature",
       y = "Frequency as top signature", x = "Sample") + 
  theme_pubr(x.text.angle = 45)

png(filename = snakemake@output[['sig_top_freq_bar_plot']], width = 12, height = 7, units = "in", res = 300)
p15 + 
  plot_annotation(title = paste0("Frequency of top ", snakemake@wildcards[['subtype']], " signatures"))
dev.off()

rm(tmp.sigtopfreq.df)



