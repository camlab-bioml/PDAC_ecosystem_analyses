suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(stringr)
})

filepaths <- snakemake@input[['sces']]
cohorts <- str_split(filepaths, "-sce-", simplify = T)[,2]
cohorts <- str_split(cohorts, ".rds", simplify = T)[,1]

dim.red <- snakemake@params[['dimred_to_plot']] %>% toupper()

print("Cohorts to be visualized together: ")
print(cohorts)

# load SCEs
sces <- lapply(filepaths, readRDS)
names(sces) <- cohorts

# subset to common genes across cohorts
sces <- lapply(sces, function(sce) {
  rownames(sce) <- paste(rowData(sce)$ensembl_id, rownames(sce), sep = "_")
  sce
})

genes <- lapply(sces, rownames)
genes.common <- Reduce(intersect, genes)
sces <- lapply(sces, function(sce) sce[genes.common,])

# tidyup metadata
sces <- lapply(sces, function(sce) {
  rowData(sce)$chr <- NULL
  rowData(sce)$gene_start <- NULL
  rowData(sce)$gene_end <- NULL
  rowData(sce)$gene_strand <- NULL
  sce
})

# cbind sces
sce <- Reduce(cbind, sces)
rm(sces)

# plot UMAPs
sce$celltype <- colData(sce)[[paste0('predicted.', snakemake@params[['annot_level']])]]
sce$celltype <- sce$celltype %>% str_to_title()

# sce$celltype <- plyr::mapvalues(sce$celltype,
#                                 from = c("Alpha", "Beta", "Delta", "Epsilon", "Gamma"),
#                                 to = c(rep("Endocrine", 5)))

p.celltype <- plotReducedDim(object = sce,
                             dimred = dim.red,
                             colour_by = 'celltype',
                             point_alpha = 0.4,
                             point_size = 0.4,
                             theme_size = 14) + 
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 0.8)))

p.cohort <- plotReducedDim(object = sce,
                           dimred = dim.red,
                           colour_by = 'cohort',
                           point_alpha = 0.4,
                           point_size = 0.4,
                           theme_size = 14) + 
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 0.8)))

# save plots
ggsave(filename = snakemake@output[['umap_celltype']], plot = p.celltype, 
       device = "png", units = "in", width = 10, height = 7, dpi = "retina", bg = "white")
ggsave(filename = snakemake@output[['umap_cohort']], plot = p.cohort, 
       device = "png", units = "in", width = 10, height = 7, dpi = "retina", bg = "white")



