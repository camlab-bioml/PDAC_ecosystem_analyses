suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(ggpubr)
  library(tidyr)
  library(stringr)
})

w <- read_tsv(snakemake@input[['gene_loading_mtx']])
#corr.thres = snakemake@params[['corr_thres']] %>% as.numeric()
#dist.thres = snakemake@params[['dist_thres']] %>% as.numeric()
celltype = snakemake@wildcards[['subtype']]

num.sim = snakemake@params[['num_sim']] %>% as.numeric()

genes <- w$gene
w$gene <- NULL
rownames(w) <- genes
rm(genes)

# get signature loadings
w.dis <- w %>% select(contains("discovery")) %>% set_rownames(., rownames(w))
w.val <- w %>% select(contains("validation")) %>% set_rownames(., rownames(w))

# simulate distance between signatures
w.dist.sim.list <- lapply(seq(length(w.dis)), function(sig.dis.id) {
  the.list <- lapply(seq(length(w.val)), function(sig.val.id) {
    sapply(seq(num.sim), function(sim.id) {
      dist(rbind(w.dis[[paste("discovery", sig.dis.id, sep = " ")]], 
                 sample(w.val[[paste("validation", sig.val.id, sep = " ")]], 
                        size = length(w.val[[paste("validation", sig.val.id, sep = " ")]]), 
                        replace = F)), 
           method = snakemake@params[['dist_method']], p = snakemake@params[['minkowski_p']] %>% as.numeric())[1]
    })
  })
  names(the.list) <- paste("validation", seq(length(w.val)), sep = " ")
  the.list
})
names(w.dist.sim.list) <- paste("discovery", seq(length(w.dis)), sep = " ")

# simulate correlation between signatures
w.corr.sim.list <- lapply(seq(length(w.dis)), function(sig.dis.id) {
  the.list <- lapply(seq(length(w.val)), function(sig.val.id) {
    sapply(seq(num.sim), function(sim.id) {
      cor(log1p(w.dis[[paste("discovery", sig.dis.id, sep = " ")]]), 
          log1p(sample(w.val[[paste("validation", sig.val.id, sep = " ")]], 
                       size = length(w.val[[paste("validation", sig.val.id, sep = " ")]]), 
                       replace = F)), 
          use = "complete.obs", method = snakemake@params[['corr_method']])
    })
  })
  names(the.list) <- paste("validation", seq(length(w.val)), sep = " ")
  the.list
})
names(w.corr.sim.list) <- paste("discovery", seq(length(w.dis)), sep = " ")
 

# save the results
saveRDS(w.dist.sim.list, file = snakemake@output[['simulated_dist_list']])
saveRDS(w.corr.sim.list, file = snakemake@output[['simulated_corr_list']])








