suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(readr)
  library(tidyr)
  library(stringr)
  library(SingleCellExperiment)
  library(NMF)
})

celltype <- snakemake@wildcards[['subtype']]
score.each.cohort <- snakemake@params[['score_each_validation_cohort']]
assay.used.for.scoring <- snakemake@params[['assay_used_for_scoring']]

cohort.correct.weight <- snakemake@params[['cohort_correct_weight']] %>% as.numeric()

# load collapsed signature gene loadings
gene_loading <- read_tsv(snakemake@input[['gene_loading_mtx']])

# load validation group gene expression
validation.exprs <- readRDS(snakemake@input[['scelist_val']])
validation.exprs <- lapply(validation.exprs, function(sce.cohort) {
  sce.cohort[intersect(gene_loading[['gene']], rownames(sce.cohort)),]
})
cohort.ids <- lapply(validation.exprs, function(sce.cohort) {
  sce.cohort$cohort
})

# run nnls
if (score.each.cohort) {
  collapsed.sig.loading.in.validation <- lapply(validation.exprs, function(sce.cohort) {
    fcnnls(x = gene_loading %>% dplyr::select(-gene) %>% as.matrix(),
           y = assay(sce.cohort, assay.used.for.scoring) %>% as.matrix())
  })
  
  collapsed.sig.loading.in.validation <- lapply(names(collapsed.sig.loading.in.validation), function(cohort) {
    cell_id <- colnames(validation.exprs[[cohort]])
    df <- cbind(collapsed.sig.loading.in.validation[[cohort]]$x %>% t() %>% as.data.frame(),
                colData(validation.exprs[[cohort]]) %>% as.data.frame() %>% select(cohort, sample))
    df$cell_id <- cell_id
    df
  })
  collapsed.sig.loading.in.validation <- Reduce(rbind, collapsed.sig.loading.in.validation)
  
  names(collapsed.sig.loading.in.validation) <- str_replace_all(names(collapsed.sig.loading.in.validation), "Rep", "RepVal")
  names(gene_loading) <- str_replace_all(names(gene_loading), "Rep", "RepVal")
  
} else{
  validation.exprs <- Reduce(cbind, validation.exprs)
  cohort.ids <- Reduce(c, cohort.ids)
  
  # construct matrix y with cohort information
  mtx.cohort.ids <- model.matrix(~ 0 + cohort_id, data = data.frame(cohort_id = cohort.ids))
  y.with.cohort.ids <- rbind(assay(validation.exprs, assay.used.for.scoring) %>% as.matrix(), 
                             mtx.cohort.ids %>% t())
  
  collapsed.sig.loading.in.validation <- fcnnls(x = rbind(gene_loading %>% dplyr::select(-gene) %>% as.matrix(),
                                                          matrix(data = cohort.correct.weight, 
                                                                 nrow = length(unique(cohort.ids)), 
                                                                 ncol = length(gene_loading)-1)),
                                                y = y.with.cohort.ids)
  
  collapsed.sig.loading.in.validation <- collapsed.sig.loading.in.validation$x %>% t() %>% as.data.frame()
  collapsed.sig.loading.in.validation <- cbind(collapsed.sig.loading.in.validation, colData(validation.exprs) %>% as.data.frame() %>% select(cohort, sample))
  collapsed.sig.loading.in.validation$cell_id <- colnames(validation.exprs)
  
  names(collapsed.sig.loading.in.validation) <- str_replace_all(names(collapsed.sig.loading.in.validation), "Rep", "RepVal")
  names(gene_loading) <- str_replace_all(names(gene_loading), "Rep", "RepVal")
  
}

# save new loading matrices
write_tsv(gene_loading, file = snakemake@output[['scored_gene_loading_mtx']])
write_tsv(collapsed.sig.loading.in.validation, file = snakemake@output[['scored_sig_loading_mtx']])
#write_tsv(signames.df, file = snakemake@output[['sig_score_namemap']])






