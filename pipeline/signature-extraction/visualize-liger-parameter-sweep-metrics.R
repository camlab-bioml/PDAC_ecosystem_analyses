suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(tidyr)
  library(readr)
  library(ggpubr)
  library(ggsci)
  library(patchwork)
})

metrics <- read_tsv(snakemake@input[['tsv']])

# calculate mean and sd across seeds
metrics <- metrics %>% 
  group_by(metric, K, Lambda) %>%
  mutate(mean = mean(overall),
         sd = sd(overall, na.rm = T)) %>%
  ungroup() %>%
  mutate(lambda = as.factor(Lambda),
         k = as.factor(K))

p1 <- ggplot(metrics, aes(x = K, y = mean, group = lambda, color = lambda)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 1, position = position_dodge(0.5)) +
  geom_line(position = position_dodge(width = 0.5)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  facet_wrap(~ metric, scales = "free_y") + 
  labs(x = "K-Value", y = "Score", title = "Selecting K") + 
  scale_color_jco() + 
  theme_pubr()

p2 <- ggplot(metrics, aes(x = Lambda, y = mean, group = k, color = k)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 1, position = position_dodge(0.05)) +
  geom_line(position = position_dodge(width = 0.05)) + 
  geom_point(position = position_dodge(width = 0.05)) +
  facet_wrap(~ metric, scales = "free_y") + 
  labs(x = "Lambda-Value", y = "Score", title = "Selecting Lambda") + 
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(metrics$k))),
                     breaks = unique(metrics$k)) + 
  theme_pubr()


ggsave(snakemake@output[["plot_k"]], p1, width = 14, height = 7, units = "in")
ggsave(snakemake@output[["plot_lambda"]], p2, width = 14, height = 7, units = "in")

