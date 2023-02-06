suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(enrichplot)
  library(patchwork)
})

ans.gse <- readRDS(snakemake@input[['gsea_go']])
ans.go <- readRDS(snakemake@input[['overrepresentation_go']])
ans.kegg <- readRDS(snakemake@input[['overrepresentation_kegg']])

condition = snakemake@wildcards[['condition']]
signature = snakemake@params[['signature']] %>% as.numeric()
if (condition == "validated") {
  w <- read_tsv(snakemake@input[['gene_loading_mtx_validated']])
  condition = snakemake@wildcards[['subtype']]
  signature = min(signature, length(grep(condition, names(w), value = T)))
} else if (condition == "collapsed") {
  w <- read_tsv(snakemake@input[['gene_loading_mtx_collapsed']])
  condition = paste(snakemake@wildcards[['subtype']], "Rep", sep = " ")
  signature = min(signature, length(grep(condition, names(w), value = T)))
}

# plot GSEA result
tab.gse <- as.data.frame(ans.gse)
tab.gse <- tab.gse %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES))

tab.gse$Description <- factor(tab.gse$Description, levels = tab.gse$Description[order(tab.gse$NES)])

p <- ggplot(tab.gse[1:20,], aes(x = NES, y = Description, fill = -log(p.adjust))) + 
  geom_bar(stat = "identity") + 
  scale_fill_viridis_c() + 
  labs(title = paste0("Signature ", signature)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))
ggsave(snakemake@output[['gsea_go_plot']], p, width = 9, height = 8, units = "in")

# plot over-representation results
p1 <- barplot(ans.go, showCategory=10) + ggtitle("GO")
ggsave(snakemake@output[['overrepresentation_go_plot']], p1, width = 7, height = 7, units = "in")

p2 <- dotplot(ans.kegg, showCategory=20) + ggtitle("KEGG")
ggsave(snakemake@output[['overrepresentation_kegg_plot']], p2, width = 7, height = 8, units = "in")






