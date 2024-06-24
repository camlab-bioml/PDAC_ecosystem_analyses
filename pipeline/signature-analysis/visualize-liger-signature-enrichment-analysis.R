suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(GSEABase)
  library(fgsea)
  library(enrichplot)
  library(patchwork)
})

# load gene ranks and 3CA pathways
geneRanks <- readRDS(snakemake@input[['gene_ranks']])
msig.C4.3ca.list <- readRDS(snakemake@input[['pathways_3ca']])

# load GSEA and ORA results
ans.gse.GO <- readRDS(snakemake@input[['gsea_go']])
ans.gse.3ca <- readRDS(snakemake@input[['gsea_3ca']])
ans.go <- readRDS(snakemake@input[['overrepresentation_go']])
ans.kegg <- readRDS(snakemake@input[['overrepresentation_kegg']])

# adjust condition and signature.id for different conditions
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
} else if (condition == "collapsed-scored-validation") {
  w <- read_tsv(snakemake@input[['gene_loading_mtx_collapsed_scored_validation']])
  condition = paste(snakemake@wildcards[['subtype']], "RepVal", sep = " ")
  signature = min(signature, length(grep(condition, names(w), value = T)))
}

# plot GSEA result
tab.gse <- as.data.frame(ans.gse.GO) %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES))

tab.gse$Description <- factor(tab.gse$Description, levels = tab.gse$Description[order(tab.gse$NES)])

p <- ggplot(tab.gse[1:20,], aes(x = NES, y = Description, fill = -log(p.adjust))) + 
  geom_bar(stat = "identity") + 
  scale_fill_viridis_c() + 
  labs(title = paste0("Signature ", signature)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))
ggsave(snakemake@output[['gsea_go_plot']], p, width = 9, height = 8, units = "in")

# plot fgsea result
topPathwaysUp <- ans.gse.3ca[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- ans.gse.3ca[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

png(snakemake@output[['gsea_3ca_plot']], width = 10, height = 7, units = "in", res = 321)
plotGseaTable(msig.C4.3ca.list[topPathways], geneRanks, ans.gse.3ca, 
              gseaParam = 0.5)
dev.off()

# plot over-representation results
p1 <- barplot(ans.go, showCategory=10) + ggtitle("GO")
ggsave(snakemake@output[['overrepresentation_go_plot']], p1, width = 7, height = 7, units = "in", dpi = 321)

tryCatch({
  p2 <- dotplot(ans.kegg, showCategory=20) + ggtitle("KEGG")
  ggsave(snakemake@output[['overrepresentation_kegg_plot']], p2, width = 7, height = 8, units = "in", dpi = 321)
}, error = function(e) {
  png(snakemake@output[['overrepresentation_kegg_plot']], width = 7, height = 8, units = "in", res = 321)
  par(mar=c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, 
       paste("doptplot() encountered this error:\n", e), 
       cex = 1, col = "black", family = "serif", font = 2, adj = 0.5)
  dev.off()
})





