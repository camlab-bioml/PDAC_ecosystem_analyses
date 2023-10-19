suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(ComplexHeatmap)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(cowplot)
})

cohort.pal <- readRDS(snakemake@input[['cohort_pal']])

w <- read_tsv(snakemake@input[['sig_top_gene_loading_mtx']])

top.genes <- read_tsv(snakemake@input[['sig_top_gene_tsv']])
sce <- readRDS(snakemake@input[['sce_subtype_dis']])
sig.loading.df <- read_tsv(snakemake@input[['sig_loading_mtx']])
plot_dir <- snakemake@params[['top_gene_exprs_plot_dir']]

ntop <- snakemake@params[["num_top_genes"]] %>% as.numeric()

w <- w |> column_to_rownames("gene")

group <- str_split(colnames(w), pattern = " Rep| RepVal", simplify = T)[,1]
row_ha <- rowAnnotation(Group = group)
w.mtx <- as.matrix(w)

png(filename = snakemake@output[['sig_top_gene_loading_heatmap']], width = ntop*ncol(w.mtx)/5, height = 7, units = "in", res = 360)
set.seed(42)
Heatmap(t(w.mtx), 
        name = "Gene loading",
        left_annotation = row_ha,
        col = viridisLite::viridis(100, option = "C"))
dev.off()

for (sig in names(top.genes)) {
  top.g <- top.genes[[sig]]
  top.g <- top.g[!is.na(top.g)]
  top.g <- top.g[top.g %in% rownames(sce)]

  print(sig)
  print(top.g)
  
  sce.top.g <- sce[top.g,]
  
  plotlist <- list()
  for (i in top.g) {
    plotlist[[i]] <- plotReducedDim(sce.top.g, dimred = "UMAP_on_Harmony", colour_by = i, by_exprs_values = "logcounts") +
       ggtitle(label = i) + 
       theme(plot.title = element_text(size = 20)) + 
       scale_color_viridis_c(option = "C")
  }

  # sig <- str_replace_all(sig, "\\.", " ")
  
  names(sig.loading.df) <- str_replace_all(names(sig.loading.df), "-|,| ", ".")

  sce.top.g$Signature_loading <- plyr::mapvalues(colnames(sce.top.g), 
                                                 from = sig.loading.df[["cell_id"]], 
                                                 to = sig.loading.df[[sig]],
                                                 warn_missing = FALSE) %>% as.numeric()

  plotlist[[sig]] <- plotReducedDim(sce.top.g, dimred = "UMAP_on_Harmony", colour_by = "Signature_loading") +
    ggtitle(label = gsub("\\.", " ", sig)) + 
    theme(plot.title = element_text(size = 20))

  plotlist[["cohort"]] <- plotReducedDim(sce.top.g, dimred = "UMAP_on_Harmony", colour_by = "cohort") +
    scale_color_manual(values = cohort.pal) +
    ggtitle(label = "Cohort") + 
    theme(plot.title = element_text(size = 20))
  
  print(paste0(sig, " plotlist length: ", length(plotlist)))

  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  png(filename = paste0(plot_dir, str_replace_all(sig, "\\.", "-"), "-top-", ntop, "-gene-exprs.png"), width = 20, height = 20, units = "in", res = 321)
  plot_grid(ncol = 3, plotlist = plotlist) %>% print()
  dev.off()
}