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

set.seed(123L)

celltypes <- snakemake@params[['celltypes']]
print(celltypes)

# load signature loading correlation dataframe
dfc.for.plot <- read_tsv(snakemake@input[["patient_profiles_correlation_data_frame"]])

print(head(dfc.for.plot))

## tidy up for plotting
dfc.for.plot <- dfc.for.plot |>
        mutate(same_cell_type_for_plot = ifelse(same_cell_type,
                "Co-occurrence within cell type",
                "Co-occurrence between cell types"
        ))

## save the data frame as tsv
write_tsv(dfc.for.plot, snakemake@output[['patient_profiles_correlation_data_frame_full_and_intra']])

# for plotting inter-cell type correlations
dfc.for.plot_list <- list()

for (ct in union(unique(dfc.for.plot$cell_type_1), unique(dfc.for.plot$cell_type_2))) {
        dfc.for.plot_list[[ct]] <- dfc.for.plot %>%
                filter(!same_cell_type) %>%
                filter(cell_type_1 == ct | cell_type_2 == ct) |>
                mutate(facet_cell_type = ct) |>
                mutate(cell_type = ifelse(cell_type_1 == ct, cell_type_2, cell_type_1))
}
rm(ct)

dfc.for.plot.inter <- Reduce(rbind, dfc.for.plot_list)
print(head(dfc.for.plot.inter))

## save the data frame as tsv
write_tsv(dfc.for.plot.inter, snakemake@output[['patient_profiles_correlation_data_frame_inter']])


# for plotting two examples of signature co-occurrance across cell types
dfc.list <- readRDS(snakemake@input[["signature_correlation_comparison_data_frame_list"]])

dfc.list <- dfc.list[c('fibroblast 4', 'CD8-positive, alpha-beta T cell 10')]

saveRDS(dfc.list, snakemake@output[['signature_correlation_comparison_data_frame_list_examples']])

# for plotting general signature co-occurrence agreement between discovery an validation
## no need to process dis.val.agree
## already done in liger-signature-correlation-comparison.R
