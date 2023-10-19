suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(readr)
  library(stringr)
})

# load validated signature gene loading correlation
gene.loading.corr.validated <- lapply(snakemake@input[['validated_gene_loading_corr']], function(path.to.file) {
  read_tsv(path.to.file)
})

celltypes <- str_split(snakemake@input[['validated_gene_loading_corr']], "signature-analysis/", simplify = T)[,2]
celltypes <- str_split(celltypes, "/signature-collapse", simplify = T)[,1]

names(gene.loading.corr.validated) <- celltypes

gene.loading.corr.validated <- lapply(gene.loading.corr.validated, function(df) {
  mtx <- df |> as.matrix()
  rownames(mtx) <- colnames(mtx)
  mtx
})

# identify validated signatures for collapsing
gene.loading.corr.validated <- lapply(gene.loading.corr.validated, function(mtx) {
  mtx > as.numeric(snakemake@params[['collapse_corr_threshold']])
})

# match up signatures for collapsing
collapse.ind <- lapply(celltypes, function(ct) {
  collapse.ind <- list(c(1))
  holder.mat <- gene.loading.corr.validated[[ct]]
  
  for (r in 1:nrow(holder.mat)) {
    for (c in 1:ncol(holder.mat)) {
      if(holder.mat[r,c]) {
        sig.num.to.find = r
        sig.num.to.add = c
        sig.num.is.added = F
        
        for (rep.sig.num in seq(length(collapse.ind))) {
          if(sig.num.to.find %in% collapse.ind[[rep.sig.num]]) {
            collapse.ind[[rep.sig.num]] <- c(collapse.ind[[rep.sig.num]], sig.num.to.add)
            sig.num.is.added = T
          }
        }
        if(!sig.num.is.added) {
          collapse.ind <- c(collapse.ind, c(sig.num.to.add))
        }
      }
    }
  }
  rm(holder.mat, r, c, sig.num.is.added, sig.num.to.add, sig.num.to.find, rep.sig.num)
  
  names(collapse.ind) <- paste(ct, "Rep", seq(length(collapse.ind)))
  
  lapply(collapse.ind, function(rep.sig) {
    unique(rep.sig)
  })
})
names(collapse.ind) <- celltypes

# make collapse guide
collapse.ind <- lapply(names(collapse.ind), function(ct) {
  ct.collapse.ind <- collapse.ind[[ct]]
  lapply(ct.collapse.ind, function(v) {
    paste(ct, v, sep = " ", collapse = " | ")
  })
})

collapse.ind <- lapply(collapse.ind, function(ct.collapse.ind) {
  data.frame(collapsed.sig.name = names(ct.collapse.ind),
             validated.sig.name = unlist(ct.collapse.ind))
})

collapse_guide <- Reduce(rbind, collapse.ind)

# save collapse guide
write_tsv(collapse_guide, file = snakemake@output[['signature_collapse_guide']])

  

