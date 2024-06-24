suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(readr)
})

cell_type_rename <- read_csv(snakemake@input[['cell_type_rename']])

# load signature interpretation
sig.interpt <- readxl::read_xlsx(snakemake@input[['sig_interpretation']])

# load cell types
celltypes <- snakemake@params[['celltypes']]
print(celltypes)

# load signature top gene analysis data frames
w.top.genes.df.list <- lapply(celltypes, function(ct) {
  w.top.genes.df.path <- grep(ct, snakemake@input[['sig_top_gene_analysis_df']], value = TRUE)[1]
  df <- read_tsv(w.top.genes.df.path)
  df$signature <- gsub(" Rep | RepVal ", " ", df$signature)
  df
})
names(w.top.genes.df.list) <- celltypes

# rbind all cell types
w.top.genes.df.combined <- Reduce(rbind, w.top.genes.df.list)

# rename signatures and cell types
w.top.genes.df.combined$signature <- plyr::mapvalues(w.top.genes.df.combined$signature,
						     from = sig.interpt$signature,
						     to = sig.interpt$`short interpretation`,
						     warn_missing = FALSE)

w.top.genes.df.combined$celltype <- plyr::mapvalues(w.top.genes.df.combined$celltype,
						    from = cell_type_rename$old_name,
						    to = cell_type_rename$new_name,
						    warn_missing = FALSE)


# remove ambient signatures
w.top.genes.df.combined <- w.top.genes.df.combined %>%
  filter(!grepl("Ambient RNA", signature))

# save
write_tsv(w.top.genes.df.combined, snakemake@output[['sig_top_gene_analysis_df_combined']])