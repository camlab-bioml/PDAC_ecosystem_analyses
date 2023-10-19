suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
})

w <- read_tsv(snakemake@input[['gene_loading_mtx']])
geneuniverse <- readRDS(snakemake@input[['geneuniverse']])

condition = snakemake@wildcards[['condition']]
signature.id = snakemake@params[['signature']] %>% as.numeric()
if (condition == "validated") {
  w <- read_tsv(snakemake@input[['gene_loading_mtx_validated']])
  condition = snakemake@wildcards[['subtype']]
  signature.id = min(signature.id, length(grep(condition, names(w), value = T)))
} else if (condition == "collapsed") {
  w <- read_tsv(snakemake@input[['gene_loading_mtx_collapsed']])
  condition = paste(snakemake@wildcards[['subtype']], "Rep", sep = " ")
  signature.id = min(signature.id, length(grep(condition, names(w), value = T)))
} else if (condition == "collapsed-scored-validation") {
  w <- read_tsv(snakemake@input[['gene_loading_mtx_collapsed_scored_validation']])
  condition = paste(snakemake@wildcards[['subtype']], "RepVal", sep = " ")
  signature.id = min(signature.id, length(grep(condition, names(w), value = T)))
}

signature <- paste(condition, signature.id, sep = " ")
print(paste("running enrichment analysis for", snakemake@wildcards[['subtype']], snakemake@wildcards[['condition']], "signature:", signature))

# get gene loadings in the signatures
genes <- w$gene
w$gene <- NULL
rownames(w) <- genes
rm(genes)

# get DE genes
degenes <- lapply(names(w), function(sig) {
  genes <- rownames(w)[w[[sig]] != 0 & !is.na(w[[sig]])]
  str_split(genes, "_", simplify = T)[,2]
})
names(degenes) <- names(w)

# GSEA using GO gene sets
geneList <- w[[signature]]
names(geneList) <- str_split(rownames(w), "_", simplify = T)[,2]
geneList <- geneList %>% sort(decreasing = T)

ans.gse <- tryCatch({
  gseGO(geneList,
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP", 
        minGSSize = 10,
        maxGSSize = as.numeric(snakemake@params[['gsea_maxGSSize']]),
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = T,
        seed = F,
        by = "fgsea")
}, error = function(e) {
  print(paste("gseGO() encountered this error:", e))
  
  gseGO(geneList,
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP", 
        minGSSize = 20,
        maxGSSize = as.numeric(snakemake@params[['gsea_maxGSSize']]),
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = T,
        seed = F,
        by = "fgsea")
})

deGenes <- degenes
geneUniverse <- geneuniverse

ans.go <- enrichGO(gene = deGenes[[signature]], 
                   ont = "BP",
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   universe = geneUniverse,
                   minGSSize = 10,
                   maxGSSize = 500,
                   readable=TRUE,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH")

# convert gene ENSEMBL ids to ENTREZ ids (may be unnecessary)
library(org.Hs.eg.db)
deGenes <- lapply(degenes, function(genes) {
  unlist(mget(genes, envir=org.Hs.egENSEMBL2EG,
              ifnotfound = NA))
})
geneUniverse <- unlist(mget(geneuniverse, envir=org.Hs.egENSEMBL2EG,
                            ifnotfound = NA))
detach("package:org.Hs.eg.db", unload = TRUE)

ans.kegg <- enrichKEGG(gene = deGenes[[signature]],
                       organism = 'hsa',
                       keyType = "kegg",
                       universe = geneUniverse,
                       pvalueCutoff = 0.05,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pAdjustMethod = "BH")

saveRDS(ans.gse, file = snakemake@output[['gsea_go']])
saveRDS(ans.go, file = snakemake@output[['overrepresentation_go']])
saveRDS(ans.kegg, file = snakemake@output[['overrepresentation_kegg']])
