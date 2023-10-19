suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
})

# compare signature correlations between collapsed and rescored validation signatures
sigprofile.corr.collapsed <- read_tsv(snakemake@input[['patient_profiles_corr_collapsed']])
sigprofile.corr.rescored.val <- read_tsv(snakemake@input[['patient_profiles_corr_collapsed_scored_val']])

sigs <- names(sigprofile.corr.collapsed)

## tidy up data
### collapsed
upper.tri.ind <- which(upper.tri(sigprofile.corr.collapsed %>% as.matrix(), diag = F), arr.ind = T)
unique.sigprofile.corr.collapsed <- cbind(upper.tri.ind, (sigprofile.corr.collapsed %>% as.matrix())[upper.tri.ind])

colnames(unique.sigprofile.corr.collapsed) <- c("sig1", "sig2", "corr")
unique.sigprofile.corr.collapsed <- unique.sigprofile.corr.collapsed %>% as.data.frame()
unique.sigprofile.corr.collapsed$sig1 <- sapply(unique.sigprofile.corr.collapsed$sig1, function(ind) names(sigprofile.corr.collapsed)[ind])
unique.sigprofile.corr.collapsed$sig2 <- sapply(unique.sigprofile.corr.collapsed$sig2, function(ind) names(sigprofile.corr.collapsed)[ind])
unique.sigprofile.corr.collapsed$condition <- "collapsed"

unique.sigprofile.corr.collapsed$sig1 <- str_replace(unique.sigprofile.corr.collapsed$sig1, "Rep ", "")
unique.sigprofile.corr.collapsed$sig2 <- str_replace(unique.sigprofile.corr.collapsed$sig2, "Rep ", "")

### rescored validation
upper.tri.ind <- which(upper.tri(sigprofile.corr.rescored.val %>% as.matrix(), diag = F), arr.ind = T)
unique.sigprofile.corr.rescored.val <- cbind(upper.tri.ind, (sigprofile.corr.rescored.val %>% as.matrix())[upper.tri.ind])

colnames(unique.sigprofile.corr.rescored.val) <- c("sig1", "sig2", "corr")
unique.sigprofile.corr.rescored.val <- unique.sigprofile.corr.rescored.val %>% as.data.frame()
unique.sigprofile.corr.rescored.val$sig1 <- sapply(unique.sigprofile.corr.rescored.val$sig1, function(ind) names(sigprofile.corr.rescored.val)[ind])
unique.sigprofile.corr.rescored.val$sig2 <- sapply(unique.sigprofile.corr.rescored.val$sig2, function(ind) names(sigprofile.corr.rescored.val)[ind])
unique.sigprofile.corr.rescored.val$condition <- "rescored validation"

unique.sigprofile.corr.rescored.val$sig1 <- str_replace(unique.sigprofile.corr.rescored.val$sig1, "RepVal ", "")
unique.sigprofile.corr.rescored.val$sig2 <- str_replace(unique.sigprofile.corr.rescored.val$sig2, "RepVal ", "")

### combine data
unique.sigprofile.corr <- data.frame(sig1 = unique.sigprofile.corr.collapsed$sig1,
                                     sig2 = unique.sigprofile.corr.collapsed$sig2,
                                     corr.collapsed = unique.sigprofile.corr.collapsed$corr,
                                     corr.rescored.val = unique.sigprofile.corr.rescored.val$corr)
unique.sigprofile.corr$celltype.pair <- paste(str_split(unique.sigprofile.corr$sig1, " ", simplify = T)[,1], 
                                              str_split(unique.sigprofile.corr$sig2, " ", simplify = T)[,1],
                                              sep = "-")

### more data manipulation
unique.sigprofile.corr <- unique.sigprofile.corr %>%
  mutate(corr.sign.collapsed = sign(corr.collapsed), corr.sign.rescored.val = sign(corr.rescored.val)) %>%
  mutate(corr.sign.equal = (corr.sign.collapsed == corr.sign.rescored.val)) %>%
  mutate(intracelltype = str_split(celltype.pair, "-", simplify = T)[,1] == str_split(celltype.pair, "-", simplify = T)[,2])

unique.sigprofile.corr.sign.count <- unique.sigprofile.corr %>%
  dplyr::count(celltype.pair, corr.sign.equal) %>%
  add_count(celltype.pair, wt = n, name = "total") %>%
  mutate(freq = n / total) %>%
  mutate(intracelltype = str_split(celltype.pair, "-", simplify = T)[,1] == str_split(celltype.pair, "-", simplify = T)[,2])

unique.sigprofile.corr.summary <- unique.sigprofile.corr %>%
  group_by(celltype.pair) %>%
  summarise(mean.corr.collapsed = mean(corr.collapsed), mean.corr.rescored.val = mean(corr.rescored.val)) %>%
  ungroup() %>%
  mutate(mean.corr.sign.collapsed = sign(mean.corr.collapsed), mean.corr.sign.rescored.val = sign(mean.corr.rescored.val)) %>%
  mutate(mean.corr.sign.equal = (mean.corr.sign.collapsed == mean.corr.sign.rescored.val)) %>%
  mutate(intracelltype = str_split(celltype.pair, "-", simplify = T)[,1] == str_split(celltype.pair, "-", simplify = T)[,2])

# save the comparison results
write_tsv(unique.sigprofile.corr %>% as.data.frame(), file = snakemake@output[['sig_loading_corr_comparison']])
write_tsv(unique.sigprofile.corr.sign.count %>% as.data.frame(), file = snakemake@output[['sig_loading_corr_sign_comparison']])
write_tsv(unique.sigprofile.corr.summary %>% as.data.frame(), file = snakemake@output[['sig_loading_corr_mean_comparison']])



# compute correlation using selected samples/cohorts ---------------------------
if(snakemake@params[['compare_per_cohort']]) {
  sigprofile.rescored.val <- read_tsv(snakemake@input[['patient_profiles_collapsed_scored_val']])
  cohorts <- sigprofile.rescored.val[['cohort']] %>% unique()
  
  individual.cohort.results.dir = snakemake@params[['comparison_results_for_individual_cohorts_dir']]
  if (!dir.exists(individual.cohort.results.dir)) dir.create(individual.cohort.results.dir)
  
  sigprofile.corr.rescored.val.list <- lapply(cohorts, function(c) {
    sigprofile.corr.rescored.val <- sigprofile.rescored.val %>%
      filter(cohort == c) %>%
      select(-cohort, -sample) %>%
      as.matrix() %>%
      cor(use = "pairwise.complete.obs") %>%
      as.data.frame()
    
    ## tidy up data
    ### rescored validation
    upper.tri.ind <- which(upper.tri(sigprofile.corr.rescored.val %>% as.matrix(), diag = F), arr.ind = T)
    unique.sigprofile.corr.rescored.val <- cbind(upper.tri.ind, (sigprofile.corr.rescored.val %>% as.matrix())[upper.tri.ind])
    
    colnames(unique.sigprofile.corr.rescored.val) <- c("sig1", "sig2", "corr")
    unique.sigprofile.corr.rescored.val <- unique.sigprofile.corr.rescored.val %>% as.data.frame()
    unique.sigprofile.corr.rescored.val$sig1 <- sapply(unique.sigprofile.corr.rescored.val$sig1, function(ind) names(sigprofile.corr.rescored.val)[ind])
    unique.sigprofile.corr.rescored.val$sig2 <- sapply(unique.sigprofile.corr.rescored.val$sig2, function(ind) names(sigprofile.corr.rescored.val)[ind])
    unique.sigprofile.corr.rescored.val$condition <- "rescored validation"
    
    unique.sigprofile.corr.rescored.val$sig1 <- str_replace(unique.sigprofile.corr.rescored.val$sig1, "RepVal ", "")
    unique.sigprofile.corr.rescored.val$sig2 <- str_replace(unique.sigprofile.corr.rescored.val$sig2, "RepVal ", "")
    
    ### combine data
    unique.sigprofile.corr <- data.frame(sig1 = unique.sigprofile.corr.collapsed$sig1,
                                         sig2 = unique.sigprofile.corr.collapsed$sig2,
                                         corr.collapsed = unique.sigprofile.corr.collapsed$corr,
                                         corr.rescored.val = unique.sigprofile.corr.rescored.val$corr)
    unique.sigprofile.corr$celltype.pair <- paste(str_split(unique.sigprofile.corr$sig1, " ", simplify = T)[,1], 
                                                  str_split(unique.sigprofile.corr$sig2, " ", simplify = T)[,1],
                                                  sep = "-")
    
    ### more data manipulation
    unique.sigprofile.corr <- unique.sigprofile.corr %>%
      mutate(corr.sign.collapsed = sign(corr.collapsed), corr.sign.rescored.val = sign(corr.rescored.val)) %>%
      mutate(corr.sign.equal = (corr.sign.collapsed == corr.sign.rescored.val)) %>%
      mutate(intracelltype = str_split(celltype.pair, "-", simplify = T)[,1] == str_split(celltype.pair, "-", simplify = T)[,2])
    
    unique.sigprofile.corr.sign.count <- unique.sigprofile.corr %>%
      dplyr::count(celltype.pair, corr.sign.equal) %>%
      add_count(celltype.pair, wt = n, name = "total") %>%
      mutate(freq = n / total) %>%
      mutate(intracelltype = str_split(celltype.pair, "-", simplify = T)[,1] == str_split(celltype.pair, "-", simplify = T)[,2])
    
    unique.sigprofile.corr.summary <- unique.sigprofile.corr %>%
      group_by(celltype.pair) %>%
      summarise(mean.corr.collapsed = mean(corr.collapsed), mean.corr.rescored.val = mean(corr.rescored.val)) %>%
      ungroup() %>%
      mutate(mean.corr.sign.collapsed = sign(mean.corr.collapsed), mean.corr.sign.rescored.val = sign(mean.corr.rescored.val)) %>%
      mutate(mean.corr.sign.equal = (mean.corr.sign.collapsed == mean.corr.sign.rescored.val)) %>%
      mutate(intracelltype = str_split(celltype.pair, "-", simplify = T)[,1] == str_split(celltype.pair, "-", simplify = T)[,2])
    
    # save the comparison results
    write_tsv(unique.sigprofile.corr, file = paste0(individual.cohort.results.dir, 
                                                    str_match(snakemake@output[['sig_loading_corr_comparison']], 
                                                              paste0("/", snakemake@wildcards[['profile']], "/\\s*(.*?)\\s*.tsv"))[,2],
                                                    "-", c, ".tsv"))
    write_tsv(unique.sigprofile.corr.sign.count, file = paste0(individual.cohort.results.dir, 
                                                               str_match(snakemake@output[['sig_loading_corr_sign_comparison']], 
                                                                         paste0("/", snakemake@wildcards[['profile']], "/\\s*(.*?)\\s*.tsv"))[,2],
                                                               "-", c, ".tsv"))
    write_tsv(unique.sigprofile.corr.summary, file = paste0(individual.cohort.results.dir, 
                                                            str_match(snakemake@output[['sig_loading_corr_mean_comparison']], 
                                                                      paste0("/", snakemake@wildcards[['profile']], "/\\s*(.*?)\\s*.tsv"))[,2],
                                                            "-", c, ".tsv"))
  })
}

