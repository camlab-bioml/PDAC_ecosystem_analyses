suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
  library(corrplot)
  library(circlize)
  library(dendextend)
})

# load the comparison results
unique.sigprofile.corr <- read_tsv(snakemake@input[['sig_loading_corr_comparison']])
unique.sigprofile.corr.sign.count <- read_tsv(snakemake@input[['sig_loading_corr_sign_comparison']])
unique.sigprofile.corr.summary <- read_tsv(snakemake@input[['sig_loading_corr_mean_comparison']])

# draw plots
ggplot(unique.sigprofile.corr %>% dplyr::filter(!intracelltype), aes(x = corr.collapsed, y = corr.rescored.val)) +
  geom_point(aes(color = celltype.pair)) +
  ggpmisc::stat_poly_line(method = "lm") +
  ggpmisc::stat_poly_eq(eq.with.lhs = "italic(hat(y))~`=`~",
                        aes(label = paste(after_stat(eq.label),
                                          after_stat(rr.label), sep = "*\", \"*"))) +
  geom_abline(color = "red", linetype = "dashed", alpha = 1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", alpha = 1) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", alpha = 1) +
  labs(title = "Each dot is correlation of a signature pair (e.g. B 1-Mono 1, Mono 3-ductal 8, CD4 2-stellate 5)") +
  theme_pubr()
ggsave(snakemake@output[['scatter_plot_all_intercell_sig_pairs']], device = "png", width = 10, height = 10, units = "in", dpi = "retina")

ggplot(unique.sigprofile.corr.summary %>% dplyr::filter(!intracelltype), aes(x = mean.corr.collapsed, y = mean.corr.rescored.val, label = celltype.pair)) +
  geom_point() +
  ggpmisc::stat_poly_line(method = "lm") +
  ggpmisc::stat_poly_eq(eq.with.lhs = "italic(hat(y))~`=`~",
                        aes(label = paste(after_stat(eq.label),
                                          after_stat(rr.label), sep = "*\", \"*"))) +
  geom_abline(color = "red", linetype = "dashed", alpha = 1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", alpha = 1) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", alpha = 1) +
  ggrepel::geom_label_repel() +
  labs(title = "Each dot is mean correlation of a cell type pair (e.g. B-Mono, Mono-ductal, CD4-stellate)") +
  theme_pubr()
ggsave(snakemake@output[['scatter_plot_intercell_sig_pair_means']], device = "png", width = 10, height = 10, units = "in", dpi = "retina")

ggplot(unique.sigprofile.corr.sign.count, aes(x = reorder(celltype.pair, -n), y = n)) +
  geom_bar(aes(fill = corr.sign.equal), stat = "identity", position = "dodge") +
  labs(title = "Whether signs of correlation between pairs of sigs are the same in collapsed and rescored validation") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(snakemake@output[['bar_plot_all_sig_pair_corr_sign_comparison']], device = "png", width = 10, height = 8, units = "in", dpi = "retina")

ggplot(unique.sigprofile.corr.sign.count %>% dplyr::filter(corr.sign.equal == TRUE), aes(x = reorder(celltype.pair, -freq), y = freq)) +
  geom_bar(aes(fill = intracelltype), stat = "identity", position = "dodge") +
  labs(title = "Freq. that pairs of sigs have the same sign for corr. between collapsed and rescored validation") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(snakemake@output[['bar_plot_sig_pair_corr_same_sign_freq']], device = "png", width = 10, height = 8, units = "in", dpi = "retina")




# compute correlation using selected samples/cohorts ---------------------------
if(snakemake@params[['compare_per_cohort']]) {
  individual.cohort.results.dir = snakemake@params[['comparison_results_for_individual_cohorts_dir']]
  individual.cohort.results.plot.dir = snakemake@params[['comparison_plots_for_individual_cohorts_dir']]
  if (!dir.exists(individual.cohort.results.plot.dir)) dir.create(individual.cohort.results.plot.dir)
  
  result.file.names <- list.files(individual.cohort.results.dir)
  
  if((length(result.file.names) > 0) & (length(result.file.names) %% 3 == 0)) {
    cohorts <- str_match(result.file.names, "comparison-\\s*(.*?)\\s*.tsv")[,2] %>% unique()
    
    for (cohort in cohorts) {
      # load the comparison results
      unique.sigprofile.corr <- read_tsv(paste0(individual.cohort.results.dir, result.file.names[grep(paste0("correlation-comparison-", cohort), result.file.names)]))
      unique.sigprofile.corr.sign.count <- read_tsv(paste0(individual.cohort.results.dir, result.file.names[grep(paste0("correlation-sign-comparison-", cohort), result.file.names)]))
      unique.sigprofile.corr.summary <- read_tsv(paste0(individual.cohort.results.dir, result.file.names[grep(paste0("correlation-mean-comparison-", cohort), result.file.names)]))
      
      # draw plots
      p1 <- ggplot(unique.sigprofile.corr %>% dplyr::filter(!intracelltype), aes(x = corr.collapsed, y = corr.rescored.val)) +
        geom_point(aes(color = celltype.pair)) +
        ggpmisc::stat_poly_line(method = "lm") +
        ggpmisc::stat_poly_eq(eq.with.lhs = "italic(hat(y))~`=`~",
                              aes(label = paste(after_stat(eq.label),
                                                after_stat(rr.label), sep = "*\", \"*"))) +
        geom_abline(color = "red", linetype = "dashed", alpha = 1) +
        geom_hline(yintercept = 0, color = "grey", linetype = "dashed", alpha = 1) +
        geom_vline(xintercept = 0, color = "grey", linetype = "dashed", alpha = 1) +
        labs(title = "Each dot is correlation of a signature pair (e.g. B 1-Mono 1, Mono 3-ductal 8, CD4 2-stellate 5)") +
        theme_pubr()
      ggsave(filename = paste0(individual.cohort.results.plot.dir, 
                               str_match(snakemake@output[['scatter_plot_all_intercell_sig_pairs']], 
                                         paste0("/", snakemake@wildcards[['profile']], "/\\s*(.*?)\\s*.png"))[,2],
                               "-", cohort, ".png"), 
             plot = p1, device = "png", width = 10, height = 10, units = "in", dpi = "retina")
      
      p2 <- ggplot(unique.sigprofile.corr.summary %>% dplyr::filter(!intracelltype), aes(x = mean.corr.collapsed, y = mean.corr.rescored.val, label = celltype.pair)) +
        geom_point() +
        ggpmisc::stat_poly_line(method = "lm") +
        ggpmisc::stat_poly_eq(eq.with.lhs = "italic(hat(y))~`=`~",
                              aes(label = paste(after_stat(eq.label),
                                                after_stat(rr.label), sep = "*\", \"*"))) +
        geom_abline(color = "red", linetype = "dashed", alpha = 1) +
        geom_hline(yintercept = 0, color = "grey", linetype = "dashed", alpha = 1) +
        geom_vline(xintercept = 0, color = "grey", linetype = "dashed", alpha = 1) +
        ggrepel::geom_label_repel() +
        labs(title = "Each dot is mean correlation of a cell type pair (e.g. B-Mono, Mono-ductal, CD4-stellate)") +
        theme_pubr()
      ggsave(filename = paste0(individual.cohort.results.plot.dir, 
                               str_match(snakemake@output[['scatter_plot_intercell_sig_pair_means']], 
                                         paste0("/", snakemake@wildcards[['profile']], "/\\s*(.*?)\\s*.png"))[,2],
                               "-", cohort, ".png"), 
             plot = p2, device = "png", width = 10, height = 10, units = "in", dpi = "retina")
      
      p3 <- ggplot(unique.sigprofile.corr.sign.count, aes(x = reorder(celltype.pair, -n), y = n)) +
        geom_bar(aes(fill = corr.sign.equal), stat = "identity", position = "dodge") +
        labs(title = "Whether signs of correlation between pairs of sigs are the same in collapsed and rescored validation") +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(filename = paste0(individual.cohort.results.plot.dir, 
                               str_match(snakemake@output[['bar_plot_all_sig_pair_corr_sign_comparison']], 
                                         paste0("/", snakemake@wildcards[['profile']], "/\\s*(.*?)\\s*.png"))[,2],
                               "-", cohort, ".png"), 
             plot = p3, device = "png", width = 10, height = 8, units = "in", dpi = "retina")
      
      p4 <- ggplot(unique.sigprofile.corr.sign.count %>% dplyr::filter(corr.sign.equal == TRUE), aes(x = reorder(celltype.pair, -freq), y = freq)) +
        geom_bar(aes(fill = intracelltype), stat = "identity", position = "dodge") +
        labs(title = "Freq. that pairs of sigs have the same sign for corr. between collapsed and rescored validation") +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(filename = paste0(individual.cohort.results.plot.dir, 
                               str_match(snakemake@output[['bar_plot_sig_pair_corr_same_sign_freq']], 
                                         paste0("/", snakemake@wildcards[['profile']], "/\\s*(.*?)\\s*.png"))[,2],
                               "-", cohort, ".png"), 
             plot = p4, device = "png", width = 10, height = 8, units = "in", dpi = "retina")
    }
  }
}

