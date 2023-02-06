suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(readr)
  library(tidyr)
  library(stringr)
})

celltype <- snakemake@wildcards[['subtype']]

# load signatures
gene_loading <- read_tsv(snakemake@input[['gene_loading_mtx']])
sig_loading <- read_tsv(snakemake@input[['sig_loading_mtx']])

# load collapse guide
collapse_guide <- read_csv(snakemake@input[['signature_collapse_guide']]) %>%
  filter(grepl(celltype, validation.sig.name))

sigs_to_collapse <- collapse_guide %>%
  filter(grepl("\\|", validated.sig.name))

# check if the signatures for these cell type need to be collapsed
if (nrow(sigs_to_collapse) > 0) {
  sigs_to_collapse <- str_split(sigs_to_collapse$validated.sig.name, " \\| ", simplify = T) %>% 
    t() %>% 
    as.data.frame() %>% 
    as.list()
  
  # extract signatures to collapse
  sigs_to_collapse.df <- lapply(sigs_to_collapse, function(l) {
    sigs <- str_subset(l, ".+")
    list(
      w = gene_loading %>% select(all_of(sigs)),
      h = sig_loading %>% select(all_of(sigs))
    )
  })
  names(sigs_to_collapse.df) <- sapply(sigs_to_collapse, function(l) l[1])
  
  # get y (explained expression)
  y.list <- lapply(sigs_to_collapse.df, function(sigs) {
    h = sigs$h %>% as.matrix()
    w = sigs$w %>% as.matrix() %>% t()
    y = h %*% w
    y
  })
  
  # get x (mean gene loading)
  x.list <- lapply(sigs_to_collapse.df, function(sigs) {
    w = sigs$w %>% as.matrix()
    x = rowMeans(w)
    x
  })
  
  # get g (collapsed signature loading)
  g.list <- lapply(seq(length(y.list)), function(sig.id) {
    y.mtx = y.list[[sig.id]]
    x = x.list[[sig.id]]
    
    apply(y.mtx, 1, function(y) {
      the.fit <- lm(y ~ 0 + x)
      coefficients(the.fit) %>% as.numeric()
    })
  })
  names(g.list) <- names(y.list)
  
  # remove old signatures
  sigs_to_remove <- sigs_to_collapse %>% unlist() %>% unname() %>% str_subset(".+")
  
  sig_loading_new <- sig_loading %>% 
    select(contains(celltype)) %>%
    select(-all_of(sigs_to_remove))
  
  gene_loading_new <- gene_loading %>%
    select(contains(celltype)) %>%
    select(-all_of(sigs_to_remove))
  
  # add collapsed signatures
  newsig_names <- collapse_guide %>%
    filter(grepl("\\|", validated.sig.name))
  newsig_names <- gsub(" ", "", newsig_names$validated.sig.name)
  
  names(g.list) <- newsig_names
  names(x.list) <- newsig_names
  
  sig_loading_new <- cbind(sig_loading_new, as.data.frame(g.list))
  gene_loading_new <- cbind(gene_loading_new, as.data.frame(x.list))
  
  # map old signature names to new signature names
  oldsig_names <- names(gene_loading_new)
  oldsig_names
  
  newsig_names <- paste(celltype, "Rep", seq(length(oldsig_names)), sep = " ")
  newsig_names
  
  signames.df <- data.frame(before.collapse = oldsig_names,
                            after.collapse = newsig_names)
  
  names(gene_loading_new) <- signames.df$after.collapse
  names(sig_loading_new) <- signames.df$after.collapse
  
  # add information to the new loading matrices
  gene_loading_new <- cbind(gene_loading_new, gene_loading %>% select(-contains(celltype)))
  sig_loading_new <- cbind(sig_loading_new, sig_loading %>% select(-contains(celltype)))
  
  # save new loading matrices and the signature name map
  write_tsv(gene_loading_new, file = snakemake@output[['collapsed_gene_loading_mtx']])
  write_tsv(sig_loading_new, file = snakemake@output[['collapsed_sig_loading_mtx']])
  write_tsv(signames.df, file = snakemake@output[['sig_collapse_namemap']])
} else {
  
  sig_loading_new <- sig_loading %>% select(contains(celltype)) 
  gene_loading_new <- gene_loading %>% select(contains(celltype)) 
  
  # map old signature names to new signature names
  oldsig_names <- names(gene_loading_new)
  oldsig_names
  
  newsig_names <- paste(celltype, "Rep", seq(length(oldsig_names)), sep = " ")
  newsig_names
  
  signames.df <- data.frame(before.collapse = oldsig_names,
                            after.collapse = newsig_names)
  
  names(gene_loading_new) <- signames.df$after.collapse
  names(sig_loading_new) <- signames.df$after.collapse
  
  # add information to the new loading matrices
  gene_loading_new <- cbind(gene_loading_new, gene_loading %>% select(-contains(celltype)))
  sig_loading_new <- cbind(sig_loading_new, sig_loading %>% select(-contains(celltype)))
  
  # save new loading matrices and the signature name map
  write_tsv(gene_loading_new, file = snakemake@output[['collapsed_gene_loading_mtx']])
  write_tsv(sig_loading_new, file = snakemake@output[['collapsed_sig_loading_mtx']])
  write_tsv(signames.df, file = snakemake@output[['sig_collapse_namemap']])
}





