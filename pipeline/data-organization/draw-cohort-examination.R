suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(stringr)
  library(sjstats)
  library(scater)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(patchwork)
})

set.seed(123L)

# load sce
sce <- readRDS(snakemake@input[['sce']])

# group cohorts
groups <- list(
  discovery = snakemake@params[["cohorts_discovery"]],
  validation = snakemake@params[["cohorts_validation"]]
)
cohorts <- unlist(groups)

# plot params
rle_style <- snakemake@params[["rle_style"]]



# set colors for cohorts
cohort.pal <- pal_npg("nrc")(length(cohorts))
names(cohort.pal) <- cohorts

# color palatte function
scale_color_cohort <- function(cohorts) {
  cohort.pal <- pal_npg("nrc")(length(cohorts))
  names(cohort.pal) <- str_split(cohorts, " ", simplify = T)[,1]
  scale_color_manual(values = cohort.pal)
}


# RLE plot for counts
p <- plotRLE(sce, exprs_values = "counts", exprs_logged = FALSE, style = rle_style, colour_by = "cohort") +
  scale_color_cohort(cohorts)

ggsave(filename = snakemake@output[['rle_plot_counts']], plot = p,
       width = 7, height = 7,
       units = "in", dpi = "retina")

# RLE plot for logcounts
p <- plotRLE(sce, exprs_values = "logcounts", style = rle_style, colour_by = "cohort") +
  scale_color_cohort(cohorts)

ggsave(filename = snakemake@output[['rle_plot_logcounts']], plot = p,
       width = 7, height = 7,
       units = "in", dpi = "retina")

print("logcounts RLE plot successfully created: ")
print(snakemake@output[["rle_plot_logcounts"]])

# RLE plot for seuratNormData
p <- plotRLE(sce, exprs_values = "seuratNormData", style = rle_style, colour_by = "cohort") +
  scale_color_cohort(cohorts)

ggsave(filename = snakemake@output[['rle_plot_seuratNormData']], plot = p,
       width = 7, height = 7,
       units = "in", dpi = "retina")

print("seuratNormData RLE plot successfully created: ")
print(snakemake@output[["rle_plot_seuratNormData"]])

# plot mean and 95% CI of logcounts coloured by cohorts
logcounts.summary.for.plot <- data.frame(cohort = sce$cohort,
                                         cell_id = colnames(sce),
                                         logcounts_mean = colMeans(logcounts(sce)),
                                         logcounts_95_perct = apply(logcounts(sce), 2, quantile, probs = 0.95, na.rm = TRUE ))

p <- ggplot(logcounts.summary.for.plot, aes(x = cell_id, y = logcounts_mean, colour = cohort)) +
  geom_point() +
  scale_color_cohort(cohorts)

ggsave(filename = snakemake@output[['logcounts_mean_plot']], plot = p,
       width = 7, height = 7,
       units = "in", dpi = "retina")

p <- ggplot(logcounts.summary.for.plot, aes(x = cell_id, y = logcounts_95_perct, colour = cohort)) +
  geom_point() +
  scale_color_cohort(cohorts)

ggsave(filename = snakemake@output[['logcounts_95_perct_plot']], plot = p,
       width = 7, height = 7,
       units = "in", dpi = "retina")