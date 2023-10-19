suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(stringr)
  library(sjstats)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(patchwork)
  library(dittoSeq)
})

set.seed(123L)

groups <- list(
  discovery = snakemake@params[["cohorts_discovery"]],
  validation = snakemake@params[["cohorts_validation"]]
)

coldata.list <- list(
  Discovery = readRDS(snakemake@input[['metadata_dis']]),
  Validation = readRDS(snakemake@input[['metadata_val']])
)

print(head(coldata.list$Discovery))
print(head(coldata.list$Validation))

cohorts <- lapply(coldata.list, function(coldata) {
  coldata$cohort %>% unique()
}) %>% unlist()

sce.list <- list(
  Discovery = readRDS(snakemake@input[["sce_dis"]]),
  Validation = readRDS(snakemake@input[["sce_val"]])
)

df.redim.list <- list(
  Discovery = readRDS(snakemake@input[['dimred_dis']]),
  Validation = readRDS(snakemake@input[['dimred_val']])
)

print(head(df.redim.list$Discovery))
print(head(df.redim.list$Validation))

# rename cell types for better presentation
cell_type_rename <- read_csv(snakemake@input[['cell_type_rename']])

df.redim.list <- lapply(df.redim.list, function(df.redim) {
  for(i in seq_len(nrow(cell_type_rename))) {
    df.redim$Cell_type <- str_replace(df.redim$Cell_type,
                                      pattern = cell_type_rename$old_name[i],
                                      replacement = cell_type_rename$new_name[i])
  }
  df.redim
})

print(head(df.redim.list$Discovery))
print(head(df.redim.list$Validation))

# set colors for cell types
celltype.pal.dis <- RColorBrewer::brewer.pal(length(unique(df.redim.list$Discovery$Cell_type)), "Set1")
names(celltype.pal.dis) <- unique(df.redim.list$Discovery$Cell_type)

celltype.pal.val <- RColorBrewer::brewer.pal(length(unique(df.redim.list$Validation$Cell_type)), "Set1")
names(celltype.pal.val) <- unique(df.redim.list$Validation$Cell_type)

celltype.pal <- data.frame(
  Cell_type_dis = names(celltype.pal.dis),
  color_dis = celltype.pal.dis,
  Cell_type_val = names(celltype.pal.val),
  color_val = celltype.pal.val
)

saveRDS(celltype.pal, snakemake@output[['celltype_pal']])

# set colors for cohorts
cohort.pal <- pal_npg("nrc")(length(cohorts))
names(cohort.pal) <- cohorts

saveRDS(cohort.pal, snakemake@output[['cohort_pal']])

# color palatte function
scale_color_cohort <- function(cohorts) {
  cohort.pal <- pal_npg("nrc")(length(cohorts))
  names(cohort.pal) <- str_split(cohorts, " ", simplify = T)[,1]
  scale_color_manual(values = cohort.pal)
}

# make metadata plots
p.list <- lapply(names(coldata.list), function(group) {
  coldata <- coldata.list[[group]]
  
  p.num_donor <- ggplot(coldata, aes(x = cohort, fill = cohort)) +
    geom_bar(stat = "count") +
    scale_fill_manual(values = cohort.pal) +
    coord_flip() +
    labs(y = "Number of samples", x = group) +
    theme_pubr() +
    theme(axis.text.y = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.title.x = element_text(face = "bold"))
  
  p.num_cell <- ggplot(coldata, aes(x = cohort, y = num_cell, color = cohort, fill = cohort)) +
    geom_violin(trim = F) +
    scale_fill_manual(values = cohort.pal) +
    scale_color_cohort(cohorts) +
    geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
    scale_y_continuous(trans = 'log10') +
    coord_flip() +
    labs(y = "Cells per sample", x = NULL) + 
    theme_pubr() +
    theme(axis.title.x = element_text(face = "bold"), axis.text.y = element_blank())
  
  p.num_gene <- ggplot(coldata, aes(x = cohort, y = num_genes, color = cohort, fill = cohort)) +
    geom_violin(trim = F) +
    scale_fill_manual(values = cohort.pal) +
    scale_color_cohort(cohorts) +
    geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
    coord_flip() +
    labs(y = "Detected genes", x = NULL) +
    theme_pubr() +
    theme(axis.title.x = element_text(face = "bold"), axis.text.y = element_blank())

  p.num_UMI <- ggplot(coldata, aes(x = cohort, y = num_UMI, color = cohort, fill = cohort)) +
    geom_violin(trim = F) +
    scale_fill_manual(values = cohort.pal) +
    scale_color_cohort(cohorts) +
    geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
    ylim(0, NA) +
    coord_flip() +
    labs(y = "Detected UMIs", x = NULL) +
    theme_pubr() +
    theme(axis.title.x = element_text(face = "bold"), axis.text.y = element_blank())
  
  p.percent_mito <- ggplot(coldata, aes(x = cohort, y = percent_mito, color = cohort, fill = cohort)) +
    geom_violin(trim = F) +
    scale_fill_manual(values = cohort.pal) +
    scale_color_cohort(cohorts) +
    geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
    coord_flip() +
    labs(y = "% Mitochondrial genes", x = NULL) +
    theme_pubr() +
    theme(axis.title.x = element_text(face = "bold"), axis.text.y = element_blank())

  p.percent_ribo <- ggplot(coldata, aes(x = cohort, y = percent_ribo, color = cohort, fill = cohort)) +
    geom_violin(trim = F) +
    scale_fill_manual(values = cohort.pal) +
    scale_color_cohort(cohorts) +
    geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
    coord_flip() +
    labs(y = "% Ribosomal genes", x = NULL) +
    theme_pubr() +
    theme(axis.title.x = element_text(face = "bold"), axis.text.y = element_blank())
  
  p.num_donor + p.num_cell + p.num_gene + p.num_UMI + p.percent_mito + p.percent_ribo + plot_layout(guides = "collect", nrow = 1) & theme(legend.position = "none")
})
names(p.list) <- names(coldata.list)
p <- ggarrange(p.list$Discovery, p.list$Validation, nrow = 2, heights = sapply(groups, function(group) {length(group)+2}))
ggsave(filename = snakemake@output[['figure1_c']], plot = p, 
       width = snakemake@params[['metadata_plot_width']], height = snakemake@params[['metadata_plot_height']], 
       units = "in", dpi = "retina")

print("Metadata plot successfully created: ")
print(snakemake@output[['figure1_c']])

# make UMAP plots
## adding number of cells per cohort in the legend
cohorts.with.num.cell <- c(df.redim.list$Discovery$Cohort, df.redim.list$Validation$Cohort) %>% unique()

cohorts
cohorts.with.num.cell
factor(str_split(cohorts.with.num.cell, " ", simplify = T)[,1], levels = cohorts) %>% rank()
cohorts.with.num.cell[order(factor(str_split(cohorts.with.num.cell, " ", simplify = T)[,1], levels = cohorts) %>% rank())]
cohorts.with.num.cell

cohort.pal <- pal_npg("nrc")(length(cohorts.with.num.cell))
names(cohort.pal) <- cohorts.with.num.cell[order(factor(str_split(cohorts.with.num.cell, " ", simplify = T)[,1], levels = cohorts) %>% rank())]

## DRAW THE UMAPs
### ggplot dotplot you can set shape = 21
p.redim.list <- lapply(df.redim.list, function(df.redim) {
  p.celltype <- ggplot(df.redim, aes(x = UMAP_1, y= UMAP_2, color = Cell_type)) +
    geom_point(alpha = 0.3, size = 0.1, shape = 1) + 
    scale_color_manual(values = c(celltype.pal.dis, celltype.pal.val)) +
    theme_pubr() +
    labs(x = "UMAP 1", y = "UMAP 2", colour = "Cell type") + 
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 0.8, shape = 16))) +
    theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
          legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))
  
  p.cohort <- ggplot(df.redim, aes(x = UMAP_1, y= UMAP_2, color = Cohort)) +
    geom_point(alpha = 0.3, size = 0.1, shape = 1) + 
    scale_color_manual(values = cohort.pal) +
    theme_pubr() +
    labs(x = "UMAP 1", y = "UMAP 2") + 
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1, shape = 16))) +
    theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
          legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))
  
  p.cohort + p.celltype + plot_layout(guides = "auto", nrow = 1) & theme(legend.position = "right")
})

p.redim <- ggarrange(p.redim.list$Discovery, p.redim.list$Validation, nrow = 2)
ggsave(filename = snakemake@output[['figure1_d']], plot = p.redim, 
       width = snakemake@params[['umap_plot_width']], height = snakemake@params[['umap_plot_height']], 
       units = "in", dpi = "retina")

print("UMAP plot successfully created")

# make cell type stacked bar plots
p.stackedbars.list <- lapply(names(df.redim.list), function(group) {
  df.redim <- df.redim.list[[group]]

  p.celltype <- ggplot(df.redim, aes(x = Sample, y = n, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c(celltype.pal.dis, celltype.pal.val)) +
    labs(x = NULL) + 
    coord_flip() +
    labs(y = "Cellt type abundance", x = NULL) +
    theme_pubr() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.position = "none"
    )

  p.celltype + plot_layout(guides = "collect", nrow = 1)
})
names(p.stackedbars.list) <- names(coldata.list)
p.stackedbars <- ggarrange(p.stackedbars.list$Discovery, p.stackedbars.list$Validation, nrow = 2)
ggsave(filename = snakemake@output[['figure1_e']], plot = p.stackedbars, 
       width = snakemake@params[['stacked_bar_plot_width']], height = snakemake@params[['stacked_bar_plot_height']], 
       units = "in", dpi = "retina")

print("Stacked bar plot successfully created")

# plot Figure 1
design <- "
AAAAAA
BBBBBC
BBBBBC
BBBBBC
"

pdf(file = snakemake@output[["figure1_pdf"]], width = snakemake@params[["figure1_width"]], height = snakemake@params[["figure1_height"]])
p + p.redim + p.stackedbars +  plot_layout(design = design)
dev.off()

print("Figure 1 PDF successfully created")

png(file = snakemake@output[["figure1_png"]], width = snakemake@params[["figure1_width"]], height = snakemake@params[["figure1_height"]], units = "in", res = 321)
p + p.redim + p.stackedbars +  plot_layout(design = design)
dev.off()

print("Figure 1 PNG successfully created")



