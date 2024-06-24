suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(magrittr)
  library(stats)
  library(gdata)
  library(posterior)
  library(bayesplot)
  library(jtools)
  library(cowplot)
  library(cmdstanr)
  library(ComplexHeatmap)
  library(PerformanceAnalytics)
  library(patchwork)
  library(ggpubr)
  library(circlize)
  library(magick)
})

cell_type_rename <- read_csv(snakemake@input[["cell_type_rename"]])

#cell_type_rename$old_name <- gsub(" ", "_", cell_type_rename$old_name)
#cell_type_rename$old_name <- gsub("-", "_", cell_type_rename$old_name)
#cell_type_rename$old_name <- gsub(",", "", cell_type_rename$old_name)

print(cell_type_rename)

# load signature interpretation
sig.interpt <- readxl::read_xlsx(snakemake@input[["sig_interpretation"]])

#sig.interpt$signature <- gsub(" ", "_", sig.interpt$signature)
#sig.interpt$signature <- gsub("-", "_", sig.interpt$signature)
#sig.interpt$signature <- gsub(",", "", sig.interpt$signature)

# signatures to remove
ambient.sigs <- read_csv(snakemake@input[["ambient_sigs"]])$signature

#ambient.sigs <- gsub(" ", "_", ambient.sigs)
#ambient.sigs <- gsub("-", "_", ambient.sigs)
#ambient.sigs <- gsub(",", "", ambient.sigs)

# load cohort color palette
cohort_pal <- readRDS(snakemake@input[["cohort_pal"]])

# load celltype color palette
celltype_pal <- readRDS(snakemake@input[["celltype_pal"]])

celltype_pal$Cell_type <- str_split(celltype_pal$Cell_type_dis, " \\(", simplify = TRUE)[, 1]

celltype_pal_to_use <- celltype_pal$color_dis
names(celltype_pal_to_use) <- celltype_pal$Cell_type

# load stan data and other parameters
celltypes <- snakemake@params[["celltypes"]]
number.of.niches <- 4
nIter <- 8000

celltypes <- gsub(" ", "_", celltypes)
celltypes <- gsub("-", "_", celltypes)
celltypes <- gsub(",", "", celltypes)

# load scematic
schematic <- image_read_pdf(snakemake@input[['schematic']]) |> image_ggplot()

# panel A
# print the schematic
ggsave(filename = snakemake@output[["figure5_a"]], plot = schematic, dpi = 360)
print("Schematic plot successfully created: ")
print(snakemake@output[["figure5_a"]])

# panel B
# draw niches
niches.to.plot <- read_tsv(snakemake@input[["microenvironment_niche_factors_combined"]]) |> column_to_rownames("niche")

## remove the unknown niches
niches.to.plot <- niches.to.plot[!(rownames(niches.to.plot) %in% c("Discovery Niche  1", "Validation Niche  2")),]

niches.to.plot <- as.matrix(niches.to.plot)

group.split <- str_split(rownames(niches.to.plot), " ", simplify = TRUE)[,1]
#ht.row.labels <- gsub("Discovery Niche|Validation Niche", "Ecotype", rownames(niches.to.plot))
row.split <- plyr::mapvalues(rownames(niches.to.plot),
  from = c("Discovery Niche  3", "Discovery Niche  2", "Discovery Niche  4", "Discovery Niche  1", "Validation Niche  3", "Validation Niche  1", "Validation Niche  4", "Validation Niche  2"),
  to = c("Ecotype 2 - Basal-like", "Ecotype 1 - Classical", "Ecotype 3 - Immune act.", "Ecotype - Unknown", "Ecotype 2 - Basal-like", "Ecotype 1 - Classical", "Ecotype 3 - Immune act.", "Ecotype - Unknown")
)

column.split <- plyr::mapvalues(gsub(" [0-9]$| [0-9][0-9]$", "", colnames(niches.to.plot)),
  from = cell_type_rename$old_name,
  to = cell_type_rename$new_name
)
colnames(niches.to.plot) <- plyr::mapvalues(colnames(niches.to.plot),
  from = sig.interpt$signature,
  to = sig.interpt$`short interpretation`
)

col_fun <- colorRamp2(seq(min(niches.to.plot), max(niches.to.plot), length.out = 100), viridisLite::viridis(100, option = "C"))

ht.niche.loadings <- Heatmap(niches.to.plot,
  height = nrow(niches.to.plot) * unit(0.3, "in"),
  row_split = row.split,
  row_title_rot = 0,
  row_title_side = "right",
  row_title_gp = gpar(fontface = "bold"),
  column_split = column.split,
  column_title_gp = gpar(fontface = "bold"),
  name = "Program loading",
  top_annotation = columnAnnotation(
    Celltype = column.split,
    col = list(Celltype = celltype_pal_to_use),
    show_annotation_name = FALSE
  ),
  left_annotation = rowAnnotation(
    Group = group.split,
    col = list(Group = c("Discovery" = "#29D9E3", "Validation" = "#B03525")),
    show_annotation_name = TRUE
  ),
  column_names_rot = 90,
  show_row_names = FALSE,
  col = col_fun
)

png(snakemake@output[["figure5_b"]], width = snakemake@params[["figure5_b_width"]], height = snakemake@params[["figure5_b_height"]], units = "in", res = 360)
draw(ht.niche.loadings#, column_title = paste0("Microenvironment Niche Factors", " - ", number.of.niches, " niches, ", nIter, " iterations")
    )
dev.off()

print("Microenvironment Niche Factors identity plot successfully created")

# panel C - REMOVED
# draw correlation bar plot
# cov.i.cor.test <- read_tsv(snakemake@input[["intrinsic_covariance_matrices_correlation_data"]])

# cell_type_rename$old_name <- gsub(" ", "_", cell_type_rename$old_name)
# cell_type_rename$old_name <- gsub("-", "_", cell_type_rename$old_name)
# cell_type_rename$old_name <- gsub(",", "", cell_type_rename$old_name)

# cov.i.cor.test$cell_type <- plyr::mapvalues(cov.i.cor.test$cell_type,
#   from = cell_type_rename$old_name,
#   to = cell_type_rename$new_name
# )

# p.cov.corr.bar <- ggplot(cov.i.cor.test, aes(x = cell_type, y = corr, fill = cell_type)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = celltype_pal_to_use) +
#   geom_hline(yintercept = 0.5, linetype = "dashed") +
#   #geom_errorbar(aes(ymin = corr - 1.96 * sqrt((1 - corr^2) / num_sig), ymax = corr + 1.96 * sqrt((1 - corr^2) / num_sig)), width = 0.2) +
#   geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
#   labs(x = "Cell type", y = "Pearson correlation\ndiscovery intrinsic cov. and validation intrinsic cov.") +
#   geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), 
#             position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt) +
#   theme_pubr() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         axis.title.x = element_blank(),
#         legend.position = "none")
# ggsave(snakemake@output[["figure5_c"]], plot = p.cov.corr.bar, width = snakemake@params[["figure5_c_width"]], height = snakemake@params[["figure5_c_height"]], units = "in", dpi = 360)

# print("Intrinsic Covariance Correlation bar plot successfully created")


# panel C
# draw niche loading scatter plot
niche.loadings <- read_tsv(snakemake@input[["microenvironment_niche_factor_loadings_to_compare"]])

niche.loadings.discovery <- niche.loadings |> filter(group == "Discovery")
niche.loadings.validation <- niche.loadings |> filter(group == "Validation")

names(niche.loadings.discovery) <- plyr::mapvalues(names(niche.loadings.discovery),
  from = c("Niche  3", "Niche  2", "Niche  4", "Niche  1"),
  to = c("Ecotype - Basal-like", "Ecotype - Classical", "Ecotype - Immune act.", "Ecotype - Immune exh.")
)
names(niche.loadings.validation) <- plyr::mapvalues(names(niche.loadings.validation),
  from = c("Niche  3", "Niche  1", "Niche  4", "Niche  2"),
  to = c("Ecotype - Basal-like", "Ecotype - Classical", "Ecotype - Immune act.", "Ecotype - Immune exh.")
)

normalize = function(v) (v - min(v, na.rm = TRUE)) / diff(range(v, na.rm = TRUE))

niche.loadings.discovery <- niche.loadings.discovery |> select(sample, "Ecotype - Basal-like", "Ecotype - Classical", "Ecotype - Immune act.", "Ecotype - Immune exh.", group) |>
  mutate("Ecotype - Basal-like" = normalize(`Ecotype - Basal-like`),
         "Ecotype - Classical" = normalize(`Ecotype - Classical`),
         "Ecotype - Immune act." = normalize(`Ecotype - Immune act.`),
         "Ecotype - Immune exh." = normalize(`Ecotype - Immune exh.`)) #|>
  #filter(grepl("PDAC", sample))
niche.loadings.validation <- niche.loadings.validation |> select(sample, "Ecotype - Basal-like", "Ecotype - Classical", "Ecotype - Immune act.", "Ecotype - Immune exh.", group) |>
  mutate("Ecotype - Basal-like" = normalize(`Ecotype - Basal-like`),
         "Ecotype - Classical" = normalize(`Ecotype - Classical`),
         "Ecotype - Immune act." = normalize(`Ecotype - Immune act.`),
         "Ecotype - Immune exh." = normalize(`Ecotype - Immune exh.`))

niche.loadings.dis.val <- rbind(niche.loadings.discovery, niche.loadings.validation)

p.niche.loadings.scatter <- ggplot(niche.loadings.dis.val, aes(x = `Ecotype - Basal-like`, y = `Ecotype - Classical`, color = group)) +
  geom_point() +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(show.legend = FALSE) +
  scale_color_manual(values = c("Discovery" = "29D9E3", "Validation" = "#B03525")) +
  labs(color = "Group") +
  theme_pubr() + 
  theme(legend.position = "none")
ggsave(snakemake@output[["figure5_c"]], plot = p.niche.loadings.scatter, width = snakemake@params[["figure5_c_width"]], height = snakemake@params[["figure5_c_height"]], units = "in", dpi = 360)

print("Niche Loadings scatter plot successfully created")

# Panel D
# load panel D
panel_d <- image_read_pdf(snakemake@input[['panel_d']]) |> image_ggplot()
ggsave(snakemake@output[["figure5_d"]], plot = panel_d, width = snakemake@params[["figure5_d_width"]], height = snakemake@params[["figure5_d_height"]], units = "in", dpi = 360)

print("Niche scores in bulk RNA-seq data plot successfully created")

# Panel E
# load panel E
panel_e <- image_read_pdf(snakemake@input[['panel_e']]) |> image_ggplot()
ggsave(snakemake@output[["figure5_e"]], plot = panel_e, width = snakemake@params[["figure5_e_width"]], height = snakemake@params[["figure5_e_height"]], units = "in", dpi = 360)

print("KRAS variants vs. niche scores plot successfully created")

# panel E - REMOVED
# plot the correlation between intrinsic covariances
cov.i.ct.dis.val <- readRDS(snakemake@input[["intrinsic_covariance_matrices_list_to_plot"]])

cov.i.ct.dis.mtx <- cov.i.ct.dis.val[["Discovery"]] |> as.matrix()
cov.i.ct.val.mtx <- cov.i.ct.dis.val[["Validation"]] |> as.matrix()

colnames(cov.i.ct.dis.mtx) <- plyr::mapvalues(colnames(cov.i.ct.dis.mtx),
  from = sig.interpt$signature,
  to = sig.interpt$`short interpretation`,
  warn_missing = FALSE
)

# calculate the Frobenius norm of the difference between the two matrices
frobenius.norm <- function(m1, m2) {
  sqrt(sum((m1 - m2)^2))
}

cov.i.ct.dis.val.f.norm <- frobenius.norm(cov.i.ct.dis.mtx, cov.i.ct.val.mtx)
cov.i.ct.dis.val.mean <- (cov.i.ct.dis.mtx + cov.i.ct.val.mtx) / 2

#col.cov = colorRamp2(seq(from = min(cov.i.ct.dis.mtx), to = max(cov.i.ct.dis.mtx), length.out = 100), viridisLite::viridis(100, option = "D"))
#col.dist = colorRamp2(seq(from = min(w.dist.mtx), to = max(w.dist.mtx), length.out = 100), rev(viridisLite::viridis(100, option = "E")))
col.cov = colorRamp2(c(-1, 0, 1), c("blue", "white", "yellow"))

ht.cov.i.dis.val <- Heatmap(cov.i.ct.dis.mtx, rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE, 
        cluster_rows = FALSE, cluster_columns = FALSE,
        layer_fun = function(j, i, x, y, w, h, fill) {
                l = i > j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.cov(pindex(cov.i.ct.dis.mtx, i[l], j[l])), col = NA))
                l = i < j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.cov(pindex(cov.i.ct.val.mtx, i[l], j[l])), col = NA))
                l = i == j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.cov(pindex(cov.i.ct.dis.val.mean, i[l], j[l])), col = "black", lwd = 3))
                #grid.rect(x[l], y[l], w[l], h[l], gp = gpar(fill = "transparent", lwd = 2))
        },
        #column_names_rot = 45,
        show_row_names = FALSE,
        row_title = "Discovery",
        row_title_side = "left",
        row_title_gp = gpar(face = "bold"),
        column_title = "Validation",
        column_title_side = "top",
        column_title_gp = gpar(face = "bold"),
        height = nrow(cov.i.ct.dis.mtx) * unit(0.7, "in"),
        #width = ncol(w.corr.mtx) * unit(0.4, "in")
)

# plot Figure 5
design <- "
AAAAAAAA
AAAAAAAA
BBBBBBBB
BBBBBBBB
BBBBBBBB
CCDDDEEE
CCDDDEEE
"

# p.cov.corr.bar <- p.cov.corr.bar + theme(legend.position = "none")
ht.niche.loadings.grob <- grid.grabExpr(draw(ht.niche.loadings, merge_legend = TRUE))

ht.cov.i.dis.val.grob = grid.grabExpr(draw(ht.cov.i.dis.val, heatmap_legend_list = list(
        Legend(title = "Intrinsic\ncovariance", col_fun = col.cov)#,
        #Legend(title = "Minkowski\nDistance", col_fun = col.dist),
        #Legend(title = "Correlation - Distance", col_fun = col.corrminusdist)
)))

pdf(file = snakemake@output[["figure5_pdf"]], width = snakemake@params[["figure5_width"]], height = snakemake@params[["figure5_height"]])
schematic + ht.niche.loadings.grob + p.niche.loadings.scatter + panel_d + panel_e + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 20))
#schematic + ht.niche.loadings.grob + p.niche.loadings.scatter + p.cov.corr.bar + ht.cov.i.dis.val.grob + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))
dev.off()

print("Figure 5 PDF successfully created")

png(file = snakemake@output[["figure5_png"]], width = snakemake@params[["figure5_width"]], height = snakemake@params[["figure5_height"]], units = "in", res = 360)
schematic + ht.niche.loadings.grob + p.niche.loadings.scatter + panel_d + panel_e + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 20))
#schematic + ht.niche.loadings.grob + p.niche.loadings.scatter + p.cov.corr.bar + ht.cov.i.dis.val.grob + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))
dev.off()

print("Figure 5 PNG successfully created")