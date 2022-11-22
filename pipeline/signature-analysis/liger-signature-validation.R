suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(stringr)
  library(corrplot)
  library(patchwork)
})

w <- read_tsv(snakemake@input[['gene_loading_mtx']])
corr.thres = snakemake@params[['corr_thres']]
dist.thres = snakemake@params[['dist_thres']]
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

# get validated signatures
sig.num = length(w.corr) / 2
w.corr.mtx <- w.corr[(sig.num+1):nrow(w.corr),1:sig.num]
sig.num = length(w.dist) / 2
w.dist.mtx <- w.dist[(sig.num+1):nrow(w.dist),1:sig.num]

validated.sig.df <- data.frame(discovery = NA, 
                               validation.corr.1 = NA, validation.1.corr = NA, validated.corr.sig = NA, validation.corr.2 = NA, validation.2.corr = NA,
                               validation.dist.1 = NA, validation.1.dist = NA, validated.dist.sig = NA, validation.dist.2 = NA, validation.2.dist = NA)
validated.corr = 1
validated.dist = 1

for (s in seq(sig.num)) {
  max.corr.1 = 0
  max.corr.2 = 0
  val.corr.1 = 0
  val.corr.2 = 0
  
  min.dist.1 = Inf
  min.dist.2 = Inf
  val.dist.1 = 0
  val.dist.2 = 0
  
  for (v in seq(sig.num)) {
    val.sig.corr = w.corr.mtx[[s]][v]
    if (val.sig.corr > max.corr.1) {
      max.corr.2 = max.corr.1
      max.corr.1 = val.sig.corr
      val.corr.2 = val.corr.1
      val.corr.1 = v
    } else if (val.sig.corr > max.corr.2) {
      max.corr.2 = val.sig.corr
      val.corr.2 = v
    }
    
    val.sig.dist = w.dist.mtx[[s]][v]
    if (val.sig.dist < min.dist.1) {
      min.dist.2 = min.dist.1
      min.dist.1 = val.sig.dist
      val.dist.2 = val.dist.1
      val.dist.1 = v
    } else if (val.sig.dist < min.dist.2) {
      min.dist.2 = val.sig.dist
      val.dist.2 = v
    }
  }
  
  holder.df <- data.frame(
    discovery = paste("discovery", s, sep = " "), 
    validation.corr.1 = paste("validation", val.corr.1, sep = " "), validation.1.corr = max.corr.1, 
    validated.corr.sig = NA, 
    validation.corr.2 = paste("validation", val.corr.2, sep = " "), validation.2.corr = max.corr.2,
    validation.dist.1 = paste("validation", val.dist.1, sep = " "), validation.1.dist = min.dist.1, 
    validated.dist.sig = NA, 
    validation.dist.2 = paste("validation", val.dist.2, sep = " "), validation.2.dist = min.dist.2
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
write_tsv(w.corr, file = snakemake@output[['gene_loading_corr']])
write_tsv(w.dist, file = snakemake@output[['gene_loading_dist']])
write_tsv(validated.sig.df, file = snakemake@output[['validated_sig_df']])
