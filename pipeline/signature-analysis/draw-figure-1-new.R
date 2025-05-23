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
  library(magick)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplotify)
})

set.seed(123L)

# load scematic
schematic <- image_read_pdf(snakemake@input[['schematic']]) |> image_ggplot()

groups <- list(
  discovery = snakemake@params[["cohorts_discovery"]],
  validation = snakemake@params[["cohorts_validation"]]
)

coldata.list <- list(
  Discovery = read_tsv(snakemake@input[['metadata_dis']]),
  Validation = read_tsv(snakemake@input[['metadata_val']])
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
  Discovery = read_tsv(snakemake@input[['dimred_dis']]),
  Validation = read_tsv(snakemake@input[['dimred_val']])
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
celltype.pal.dis <- RColorBrewer::brewer.pal(length(unique(df.redim.list$Discovery$Cell_type)), snakemake@params[["cell_type_pallete"]])
names(celltype.pal.dis) <- unique(df.redim.list$Discovery$Cell_type)

celltype.pal.val <- RColorBrewer::brewer.pal(length(unique(df.redim.list$Validation$Cell_type)), snakemake@params[["cell_type_pallete"]])
names(celltype.pal.val) <- unique(df.redim.list$Validation$Cell_type)

celltype.pal <- data.frame(
  Cell_type_dis = names(celltype.pal.dis),
  color_dis = celltype.pal.dis,
  Cell_type_val = names(celltype.pal.val),
  color_val = celltype.pal.val
)

#saveRDS(celltype.pal, snakemake@output[['celltype_pal']])

# set colors for cohorts
cohort.pal <- pal_npg("nrc")(length(cohorts))
names(cohort.pal) <- cohorts

#saveRDS(cohort.pal, snakemake@output[['cohort_pal']])

# color palatte function
scale_color_cohort <- function(cohorts) {
  cohort.pal <- pal_npg("nrc")(length(cohorts))
  names(cohort.pal) <- str_split(cohorts, " ", simplify = T)[,1]
  scale_color_manual(values = cohort.pal)
}

# print the schematic
ggsave(filename = snakemake@output[["figure1_a"]], plot = schematic, dpi = 600)
print("Schematic plot successfully created: ")
print(snakemake@output[["figure1_a"]])

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
    #scale_y_continuous(labels = scales::scientific) +
    theme_pubr() +
    theme(axis.title.x = element_text(face = "bold"), axis.text.y = element_blank())

  if (group == "Discovery") {
    p.num_cell <- p.num_cell + scale_y_continuous(trans = 'log10', breaks = c(300, 1000, 5000, 30000))
  }
  
  p.num_gene <- ggplot(coldata, aes(x = cohort, y = num_genes, color = cohort, fill = cohort)) +
    geom_violin(trim = F) +
    scale_fill_manual(values = cohort.pal) +
    scale_color_cohort(cohorts) +
    geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
    coord_flip() +
    labs(y = "Detected genes", x = NULL) +
    #scale_y_continuous(labels = scales::scientific) +
    theme_pubr() +
    theme(axis.title.x = element_text(face = "bold"), axis.text.y = element_blank())

  p.num_UMI <- ggplot(coldata, aes(x = cohort, y = num_UMI, color = cohort, fill = cohort)) +
    geom_violin(trim = F) +
    scale_fill_manual(values = cohort.pal) +
    scale_color_cohort(cohorts) +
    geom_boxplot(width = 0.1, fill = "white", color = "grey40") +
    ylim(1, NA) +
    coord_flip() +
    labs(y = "Detected UMIs", x = NULL) +
    #scale_y_continuous(labels = scales::scientific) +
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
ggsave(filename = snakemake@output[['figure1_b_supp']], plot = p, 
       width = snakemake@params[['metadata_plot_width']], height = snakemake@params[['metadata_plot_height']], 
       units = "in", dpi = 600)

print("Metadata plot successfully created: ")
print(snakemake@output[['figure1_b']])

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
  #df.redim <- df.redim.list[[group]]
  p.celltype <- ggplot(df.redim, aes(x = UMAP_1, y= UMAP_2, color = Cell_type)) +
    geom_point(alpha = 0.3, size = 0.1, shape = 1) + 
    scale_color_manual(values = c(celltype.pal.dis, celltype.pal.val)) +
    theme_pubr() +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Cell type") + 
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8, shape = 16, byrow = TRUE))) +
    theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
          legend.title = element_text(face = "bold", size = 10), legend.text = element_text(face = "bold", size = 8),
          legend.spacing.y = unit(0.1, "cm"))
  
  p.cohort <- ggplot(df.redim, aes(x = UMAP_1, y= UMAP_2, color = Cohort)) +
    geom_point(alpha = 0.3, size = 0.1, shape = 1) + 
    scale_color_manual(values = cohort.pal) +
    theme_pubr() +
    labs(x = "UMAP 1", y = "UMAP 2") + 
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1, shape = 16, byrow = TRUE))) +
    theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
          legend.title = element_text(face = "bold", size = 10), legend.text = element_text(face = "bold", size = 8),
          legend.spacing.y = unit(0.1, "cm"))
  
  p.cohort + p.celltype + plot_layout(guides = "auto", nrow = 1) & theme(legend.position = "right", legend.spacing.y = unit(0.1, "cm"))
})

p.redim <- ggarrange(p.redim.list$Discovery, p.redim.list$Validation, nrow = 2)
ggsave(plot = p.redim, filename = snakemake@output[['figure1_b']], device = "png",
       width = snakemake@params[['umap_plot_width']], height = snakemake@params[['umap_plot_height']], 
       units = "in", dpi = 600)

print("UMAP plot successfully created")

# make cell type stacked bar plots
p.stackedbars.list <- lapply(names(df.redim.list), function(group) {
  df.redim <- df.redim.list[[group]]

  p.celltype <- ggplot(df.redim, aes(x = Sample, y = n, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c(celltype.pal.dis, celltype.pal.val)) +
    #labs(x = NULL) + 
    coord_flip() +
    labs(y = "Cell type abundance", x = group) +
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
p.stackedbars <- ggarrange(p.stackedbars.list$Discovery + theme(axis.title.x = element_blank()), 
			   p.stackedbars.list$Validation, 
			   nrow = 2)
ggsave(filename = snakemake@output[['figure1_c']], plot = p.stackedbars, 
       width = snakemake@params[['stacked_bar_plot_width']], height = snakemake@params[['stacked_bar_plot_height']], 
       units = "in", dpi = 600)

print("Stacked bar plot successfully created")

# make cell type marker dot plots
p.markerdots.list <- lapply(names(sce.list), function(group) {
  sce <- sce.list[[group]]
  rownames(sce) <- str_split(rownames(sce), "_", simplify = T)[,2]
  sce$celltype <- plyr::mapvalues(sce$celltype, from = cell_type_rename$old_name, cell_type_rename$new_name)

  p.markerdots <- dittoDotPlot(sce, 
                               vars = c("EPCAM", "CDH1", "KRT19", "KRT5", "KRT8", "CDH5", "COL1A1", "COL1A2", "VIM", "ACTA2", "PTPRC", "CD68", "CD14", "FCGR3A", "CD1C", "CLEC9A", "MS4A1", "CD3D", "CD3E", "CD4", "CD8A", "NKG7"), 
                               group.by = "celltype", 
                               ylab = NULL,
                               legend.color.title = "Relative\nexpression",
                               legend.size.title = "Percent\nexpression")
})
names(p.markerdots.list) <- names(sce.list)
p.markerdots <- ggarrange(p.markerdots.list$Discovery + theme(legend.position = "none"), p.markerdots.list$Validation, nrow = 1)
ggsave(filename = snakemake@output[['figure1_d_supp']], plot = p.markerdots, 
       width = snakemake@params[['marker_dot_plot_width']], height = snakemake@params[['marker_dot_plot_height']], 
       units = "in", dpi = 600)

print("Marker dot plot successfully created")

################################### NEW STUFF ###################################
# load signature interpretation
sig.interpt <- readxl::read_xlsx(snakemake@input[['sig_interpretation']])
# signatures to remove
ambient.sigs <- read_csv(snakemake@input[['ambient_sigs']])$signature

# panel D
# load signature numbers df
sig.number.df <- readRDS(snakemake@input[['sig_number']])

sig.number.df <- sig.number.df |>
  mutate(cell_type_to_show = plyr::mapvalues(cell_type,
                                             from = cell_type_rename$old_name,
                                             to = cell_type_rename$new_name))

sig.number.df.purged <- data.frame(
  Condition = "Contextualized",
  cell_type = cell_type_rename$old_name,
  Number = c(8, 3, 4, 10, 5, 5, 11, 4),
  cell_type_to_show = cell_type_rename$new_name
)

sig.number.df <- rbind(sig.number.df, sig.number.df.purged)

# plot program numbers
p.signum <- ggplot(sig.number.df, 
    aes(x = factor(Condition, levels = c("Discovery", "Validated", "Merged", "Contextualized")),
        y = 1,
        fill = Number)) +
  geom_tile(color = "white", width = 1) +
  geom_text(aes(label = Number), color = "white") +
  scale_fill_viridis_c(name = "Number of\nprograms") +
  scale_x_discrete(expand = c(0, 0)) +
  facet_wrap(~ cell_type_to_show, nrow = 1) +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = "none",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 0, -10, 0)
  )

ggsave(snakemake@output[['figure1_d']], width = 12, height = 1.5, units = "in", dpi = 600)
print("Program number plot successfully created")

# panel D
gene.loading.top.gene.df.to.plot.list <- readRDS(snakemake@input[['selected_sig_top_gene_loading_mtx_list']])

gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.to.plot.list$`pancreatic epithelial cell`
ht.epi <- Heatmap(gene.loading.top.gene.df.to.plot |> t(),
                  name = "Gene loading",
                  #width = unit(snakemake@params[['figure2_d_width']] - 3, "in"),
                  height = unit(0.95, "in"),
                  row_labels = plyr::mapvalues(gsub(" Rep ", " ", colnames(gene.loading.top.gene.df.to.plot)),
                                               from = sig.interpt$signature,
                                               to = sig.interpt$`short interpretation`,
                                               warn_missing = FALSE),
                  column_names_side = "bottom",
                  #column_names_rot = 45,
                  column_names_gp = gpar(fontsize = 8),
                  column_dend_side = "top",
		              #show_column_den = FALSE,
		              show_row_dend = FALSE,
                  col = colorRamp2(seq(from = 0, to = 0.7, length.out = 100), viridisLite::viridis(100, option = "C")),
                  show_heatmap_legend = FALSE)

gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.to.plot.list$`fibroblast`
ht.fibro <- Heatmap(gene.loading.top.gene.df.to.plot |> t(),
                    name = "Gene loading",
                    #width = unit(snakemake@params[['figure2_d_width']] - 3, "in"),
                    height = unit(0.8, "in"),
                    row_labels = plyr::mapvalues(gsub(" Rep ", " ", colnames(gene.loading.top.gene.df.to.plot)),
                                                 from = sig.interpt$signature,
                                                 to = sig.interpt$`short interpretation`,
                                                 warn_missing = FALSE),
                    column_names_side = "bottom",
                    #column_names_rot = 45,
                    column_names_gp = gpar(fontsize = 8),
                    column_dend_side = "top",
		                #show_column_dend = FALSE,
		                show_row_dend = FALSE,
                    col = colorRamp2(seq(from = 0, to = 0.7, length.out = 100), viridisLite::viridis(100, option = "C")),
                    show_heatmap_legend = FALSE)

png(snakemake@output[['figure1_e']], width = 12, height = 3, units = "in", res = 600)
ggarrange(as.grob(ht.epi), nrow = 1)
dev.off()

png(snakemake@output[['figure1_f']], width = 12, height = 3, units = "in", res = 600)
ggarrange(as.grob(ht.fibro), nrow = 1)
dev.off()

print("Selected signature top marker plot successfully created")

################################### NEW STUFF ###################################
# plot Figure D
design <- "
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
BBBBBBBCC
BBBBBBBCC
BBBBBBBCC
BBBBBBBCC
BBBBBBBCC
BBBBBBBCC
DDDDDDDDD
EEEEEEEEE
EEEEEEEEE
EEEEEEEEE
FFFFFFFFF
FFFFFFFFF
"

pdf(file = snakemake@output[["figure1_pdf"]], width = snakemake@params[["figure1_width"]], height = snakemake@params[["figure1_height"]])
schematic + p.redim + p.stackedbars + p.signum + as.grob(ht.epi) + as.grob(ht.fibro) + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()

print("Figure 1 PDF successfully created")

png(file = snakemake@output[["figure1_png"]], width = snakemake@params[["figure1_width"]], height = snakemake@params[["figure1_height"]], units = "in", res = 600)
schematic + p.redim + p.stackedbars + p.signum + as.grob(ht.epi) + as.grob(ht.fibro) + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()

print("Figure 1 PNG successfully created")



