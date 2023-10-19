suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(stringr)
  library(clue)
  library(corrplot)
  library(patchwork)
})

w <- read_tsv(snakemake@input[['gene_loading_mtx']])
corr.thres = snakemake@params[['corr_thres']] %>% as.numeric()
dist.thres = snakemake@params[['dist_thres']] %>% as.numeric()
celltype = snakemake@wildcards[['subtype']]

genes <- w$gene
w$gene <- NULL
rownames(w) <- genes
rm(genes)

# compute correlation between discovery and validation
w.corr <- cor(log1p(w), use = "complete.obs", method = snakemake@params[['corr_method']]) %>%
  as.data.frame()

# compute distance between discovery and validation
w.mtx <- as.matrix(w) %>% t()
w.dist <- dist(w.mtx, method = snakemake@params[['dist_method']], p = as.numeric(snakemake@params[['minkowski_p']])) %>% 
  as.matrix() %>%
  as.data.frame()

# Hungarian method - forcing one-to-one mapping
sig.num = length(w.corr) / 2
w.corr.mtx <- w.corr[1:sig.num, (sig.num+1):ncol(w.corr)] %>% as.matrix()
w.corr.mtx[w.corr.mtx < 0] = 0
sig.num = length(w.dist) / 2
w.dist.mtx <- w.dist[1:sig.num, (sig.num+1):ncol(w.dist)] %>% as.matrix()

map.corr <- solve_LSAP(w.corr.mtx, maximum = T)
map.dist <- solve_LSAP(w.dist.mtx, maximum = F)

w.dist.mtx <- w.dist.mtx[,map.dist]
w.corr.mtx <- w.corr.mtx[,map.corr]

## get validated signatures
validated.sig.df <- data.frame(discovery = NA, 
                               validation.corr.1 = NA, validation.1.corr = NA, validated.corr.sig = NA, validation.corr.2 = NA, validation.2.corr = NA,
                               validation.dist.1 = NA, validation.1.dist = NA, validated.dist.sig = NA, validation.dist.2 = NA, validation.2.dist = NA)
validated.corr = 1
validated.dist = 1

for (s in seq(sig.num)) {
  max.corr.1 = w.corr.mtx[s,s]
  max.corr.2 = max(w.corr.mtx[s,][-s])
  val.corr.1 = colnames(w.corr.mtx)[s]
  val.corr.2 = (colnames(w.corr.mtx)[-s])[which.max(w.corr.mtx[s,][-s])]
  
  min.dist.1 = w.dist.mtx[s,s]
  min.dist.2 = min(w.dist.mtx[s,][-s])
  val.dist.1 = colnames(w.dist.mtx)[s]
  val.dist.2 = (colnames(w.dist.mtx)[-s])[which.min(w.dist.mtx[s,][-s])]
  
  holder.df <- data.frame(
    discovery = rownames(w.corr.mtx)[s], 
    validation.corr.1 = val.corr.1, validation.1.corr = max.corr.1, 
    validated.corr.sig = NA, 
    validation.corr.2 = val.corr.2, validation.2.corr = max.corr.2,
    validation.dist.1 = val.dist.1, validation.1.dist = min.dist.1, 
    validated.dist.sig = NA, 
    validation.dist.2 = val.dist.2, validation.2.dist = min.dist.2
  )
  if (max.corr.1 > corr.thres) {
    holder.df$validated.corr.sig <- paste(celltype, validated.corr, sep = " ")
    validated.corr = validated.corr + 1
  }
  if (min.dist.1 < dist.thres) {
    holder.df$validated.dist.sig <- paste(celltype, validated.dist, sep = " ")
    validated.dist = validated.dist + 1
  }
  print(holder.df)
  validated.sig.df <- rbind(validated.sig.df, holder.df)
}

validated.sig.df <- validated.sig.df %>% filter(!if_all(everything(), is.na))

# save the results
#write_tsv(w.corr, file = snakemake@output[['gene_loading_corr']])
#write_tsv(w.dist, file = snakemake@output[['gene_loading_dist']])
write_tsv(validated.sig.df, file = snakemake@output[['validated_sig_df']])
