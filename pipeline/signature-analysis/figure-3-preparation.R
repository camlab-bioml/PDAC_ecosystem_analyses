suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(readr)
  library(BiocParallel)
  library(sjstats)
  library(dittoSeq)
  library(ggplotify)
  library(patchwork)
  library(cowplot)
  library(circlize)
  library(scales)
})

celltypes <- snakemake@params[['celltypes']]
print(celltypes)

# load signature loading variance dataframes
sig.loading.var.dis <- read_tsv(snakemake@input[['sig_loading_var_dis']]) |> column_to_rownames("signature")
sig.loading.var.val <- read_tsv(snakemake@input[['sig_loading_var_val']]) |> column_to_rownames("signature")

print(data.frame(discovery = sig.loading.var.dis$signature, validation = sig.loading.var.val$signature))

sig.split <- str_split(rownames(sig.loading.var.dis), " Sig ", simplify = TRUE)[, 1]
unique(sig.split)

names(sig.loading.var.dis) <- c(
        "Intra-patient heterogeneity Discovery",
        "Inter-patient heterogeneity Discovery (mean)",
        "Inter-patient heterogeneity Discovery (median)")
names(sig.loading.var.val) <- c(
        "Intra-patient heterogeneity Validation",
        "Inter-patient heterogeneity Validation (mean)",
        "Inter-patient heterogeneity Validation (median)")

sig.loading.var.dis.val <- full_join(sig.loading.var.dis |> rownames_to_column("signature"),
                                     sig.loading.var.val |> rownames_to_column("signature"),
                                     by = "signature")

print(head(sig.loading.var.dis.val))

write_tsv(sig.loading.var.dis.val, snakemake@output[['sig_loading_var_dis_val']])