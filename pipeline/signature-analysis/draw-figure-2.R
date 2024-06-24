suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(stringr)
  library(scales)
  library(BiocParallel)
  library(sjstats)
  library(dittoSeq)
  library(ggplotify)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(cowplot)
  library(circlize)
  library(ComplexHeatmap)
})

cohort.pal <- readRDS(snakemake@input[["cohort_pal"]])

cell_type_rename <- read_csv(snakemake@input[['cell_type_rename']])

# load signature interpretation
sig.interpt <- readxl::read_xlsx(snakemake@input[['sig_interpretation']])

# signatures to remove
ambient.sigs <- read_csv(snakemake@input[['ambient_sigs']])$signature

# panel A
# load signature numbers df
sig.number.df <- readRDS(snakemake@input[['sig_number']])

sig.number.df <- sig.number.df |>
  mutate(cell_type_to_show = plyr::mapvalues(cell_type,
                                             from = cell_type_rename$old_name,
                                             to = cell_type_rename$new_name))

# plot signature numbers
p.signum <- ggplot(sig.number.df, aes(x = factor(Condition, levels = c("Merged", "Validated", "Discovery")), y = Number, fill = Condition)) +
  geom_bar(stat = 'identity', position = "dodge") +
  scale_fill_npg() +
  facet_wrap(~cell_type_to_show, nrow = 2) +
  labs(#title = paste0('Number of signatures for each cell types'),
    x = NULL, 
    y = 'Number of signatures') +
  coord_flip() +
  ggpubr::theme_pubr() #+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(snakemake@output[['figure2_a']], width = snakemake@params[['figure2_a_width']], height = snakemake@params[['figure2_a_height']], units = "in", dpi = 360)

print("Signature number plot successfully created")

# panel B
w.corr <- read_tsv(snakemake@input[["gene_loading_corr"]])
rownames(w.corr) <- names(w.corr)

w.dist <- read_tsv(snakemake@input[["gene_loading_dist"]])
rownames(w.dist) <- names(w.dist)

sig.num <- length(w.corr) / 2
w.corr.mtx <- w.corr[1:sig.num, (sig.num + 1):ncol(w.corr)] %>% as.matrix()
w.corr.mtx[w.corr.mtx < 0] <- 0
sig.num <- length(w.dist) / 2
w.dist.mtx <- w.dist[1:sig.num, (sig.num + 1):ncol(w.dist)] %>% as.matrix()

validated.sig.df <- read_tsv(snakemake@input[["validated_sig_df"]])

map.corr <- str_split(validated.sig.df[["validation.corr.1"]], " ", simplify = T)[, 2] %>% as.numeric()
map.dist <- str_split(validated.sig.df[["validation.dist.1"]], " ", simplify = T)[, 2] %>% as.numeric()

w.corr.mtx <- w.corr.mtx[, map.corr]
w.dist.mtx <- w.dist.mtx[, map.corr]

w.corr.mtx <- matrix(as.numeric(w.corr.mtx),
        ncol = ncol(w.corr.mtx),
        dimnames = list(
                paste("Discovery", seq_len(nrow(w.corr.mtx)), sep = " "),
                str_to_title(colnames(w.corr.mtx))
        )
)
w.dist.mtx <- matrix(as.numeric(w.dist.mtx) %>% rescale(),
        ncol = ncol(w.dist.mtx),
        dimnames = list(
                paste("Discovery", seq_len(nrow(w.dist.mtx)), sep = " "),
                str_to_title(colnames(w.dist.mtx))
        )
)

w.corrminusdist.mtx <- w.corr.mtx - w.dist.mtx

col.corr = colorRamp2(seq(from = min(w.corr.mtx), to = max(w.corr.mtx), length.out = 100), viridisLite::viridis(100, option = "D"))
col.dist = colorRamp2(seq(from = min(w.dist.mtx), to = max(w.dist.mtx), length.out = 100), rev(viridisLite::viridis(100, option = "E")))
col.corrminusdist = colorRamp2(c(-1, 0, 1), c("blue", "white", "yellow"))

# plot gene loading correlations between discovery and validation
ht.coor <- Heatmap(w.corr.mtx,
        name = "Correlation", rect_gp = gpar(type = "none"),
        width = ncol(w.corr.mtx) * unit(0.4, "in"),
        cluster_rows = F, cluster_columns = F,
        col = colorRamp2(seq(from = min(w.corr.mtx), to = max(w.corr.mtx), length.out = 100), viridisLite::viridis(100, option = "D")),
        # cell_fun = function(j, i, x, y, w, h, fill) {
        #         if (i > j) {
        #                 grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "white"))
        #         }
        #         if (i == j) {
        #                 grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "black", lwd = 2))
        #         }
        # },
        layer_fun = function(j, i, x, y, w, h, fill) {
                l = i > j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.corr(pindex(w.corr.mtx, i[l], j[l])), col = NA))
                l = i == j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.corr(pindex(w.corr.mtx, i[l], j[l])), col = "black", lwd = 2))
        },
        show_row_names = FALSE
)

# plot gene loading distance between discovery and validation
ht.dist <- Heatmap(w.dist.mtx,
        name = "Distance", rect_gp = gpar(type = "none"),
        width = ncol(w.dist.mtx) * unit(0.4, "in"),
        cluster_rows = F, cluster_columns = F,
        col = colorRamp2(seq(from = min(w.dist.mtx), to = max(w.dist.mtx), length.out = 100), rev(viridisLite::viridis(100, option = "E"))),
        # cell_fun = function(j, i, x, y, w, h, fill) {
        #         if (i < j) {
        #                 grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "white"))
        #         }
        #         if (i == j) {
        #                 grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "black", lwd = 2))
        #         }
        # },
        layer_fun = function(j, i, x, y, w, h, fill) {
                l = i < j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.dist(pindex(w.dist.mtx, i[l], j[l])), col = NA))
                l = i == j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.dist(pindex(w.dist.mtx, i[l], j[l])), col = "black", lwd = 2))
        },
        show_column_names = FALSE
)

ht.dist.corr <- Heatmap(w.corr.mtx, rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE, 
        cluster_rows = FALSE, cluster_columns = FALSE,
        layer_fun = function(j, i, x, y, w, h, fill) {
                l = i > j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.corr(pindex(w.corr.mtx, i[l], j[l])), col = NA))
                l = i < j
                grid.rect(x[l], y[l], w[l], h[l],
                        gp = gpar(fill = col.dist(pindex(w.dist.mtx, i[l], j[l])), col = NA))
                l = i == j
                grid.rect(x[l], y[l], w[l], h[l],
                        #gp = gpar(fill = col.corrminusdist(pindex(w.corrminusdist.mtx, i[l], j[l])), col = "black", lwd = 3)
                        gp = gpar(fill = col.corr(pindex(w.corr.mtx, i[l], j[l])), col = "black", lwd = 3))
                #grid.rect(x[l], y[l], w[l], h[l], gp = gpar(fill = "transparent", lwd = 2))
        },
        #column_names_rot = 45,
        height = nrow(w.corr.mtx) * unit(0.3, "in"),
        width = ncol(w.corr.mtx) * unit(0.4, "in"))

# ht.dist.corr <- Heatmap(w.corr.mtx, rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE, 
#         cluster_rows = FALSE, cluster_columns = FALSE,
#         cell_fun = function(j, i, x, y, w, h, fill) {
#                 if (i > j) {
#                         grid.rect(x, y, w, h, gp = gpar(fill = col.corr(pindex(w.corr.mtx, i, j)), col = NA))
#                 } else if (i == j) {
#                         grid.rect(x, y, w, h, gp = gpar(fill = col.corrminusdist(pindex(w.corrminusdist.mtx, i, j)), col = "black", lwd = 2))
#                 } else {
#                         grid.rect(x, y, w, h, gp = gpar(fill = col.dist(pindex(w.dist.mtx, i, j)), col = NA))
#                 }
#         },
#         #column_names_rot = 45,
#         height = nrow(w.corr.mtx) * unit(0.3, "in"),
#         width = ncol(w.corr.mtx) * unit(0.4, "in"))

# plot both gene loading correlations and distance between discovery and validation
png(snakemake@output[["figure2_b"]], width = snakemake@params[["figure2_b_width"]], height = snakemake@params[["figure2_b_height"]], units = "in", res = 360)
draw(ht.coor + ht.dist, ht_gap = unit(-80, "mm"))
dev.off()

# panel C
# heatmap annotations
sig.loading.var.dis <- read_tsv(snakemake@input[['sig_loading_var_dis']]) |> column_to_rownames("signature")
sig.loading.var.val <- read_tsv(snakemake@input[['sig_loading_var_val']]) |> column_to_rownames("signature")

print(data.frame(discovery = sig.loading.var.dis$signature, validation = sig.loading.var.val$signature))

sig.loading.var.to.plot.dis <- sig.loading.var.dis[!(rownames(sig.loading.var.dis) %in% ambient.sigs),]
sig.loading.var.to.plot.val <- sig.loading.var.val[!(rownames(sig.loading.var.val) %in% ambient.sigs),]

sig.split <- str_split(rownames(sig.loading.var.to.plot.dis), " Sig ", simplify = TRUE)[, 1]
unique(sig.split)
sig.split <- plyr::mapvalues(sig.split,
        from = cell_type_rename$old_name,
        to = cell_type_rename$new_name)

names(sig.loading.var.to.plot.dis) <- c(
        "Intra-patient\nheterogeneity\nDiscovery",
        "Inter-patient\nheterogeneity\nDiscovery (mean)",
        "Inter-patient\nheterogeneity\nDiscovery")
names(sig.loading.var.to.plot.val) <- c(
        "Intra-patient\nheterogeneity\nValidation",
        "Inter-patient\nheterogeneity\nValidation (mean)",
        "Inter-patient\nheterogeneity\nValidation")

sig.loading.var.to.plot <- full_join(sig.loading.var.to.plot.dis |> rownames_to_column("signature"),
                                     sig.loading.var.to.plot.val |> rownames_to_column("signature"),
                                     by = "signature")

print(head(sig.loading.var.to.plot))

# mtrices for plotting
sig.loading.mtx.df.list <- list(
	"collapsed" = read_tsv(snakemake@input[['sig_loading_median_mtx_dis']]),
	"collapsed-scored-validation" = read_tsv(snakemake@input[['sig_loading_median_mtx_val']])
)

rescale_1_99 <- function(x) {
        (x - quantile(x, probs = c(0.01), na.rm = TRUE)) /
                (quantile(x, probs = c(0.99), na.rm = TRUE) - quantile(x, probs = c(0.01), na.rm = TRUE))
}

mtx.to.plot.dis <- sig.loading.mtx.df.list[["collapsed"]] |>
        select(-all_of(ambient.sigs)) |>
        reframe(across(where(is.numeric), rescale)) |>
        as.matrix()
mtx.to.plot.val <- sig.loading.mtx.df.list[["collapsed-scored-validation"]] |>
        select(-all_of(ambient.sigs)) |>
        reframe(across(where(is.numeric), rescale)) |>
        as.matrix()

# colnames(mtx.to.plot.dis) <- paste0(sig.split, " Sig ", str_split(colnames(mtx.to.plot.dis), " Sig ", simplify = TRUE)[,2])
# colnames(mtx.to.plot.val) <- paste0(sig.split, " Sig ", str_split(colnames(mtx.to.plot.val), " Sig ", simplify = TRUE)[,2])

colnames(mtx.to.plot.dis) <- plyr::mapvalues(gsub(" Sig ", " ", colnames(mtx.to.plot.dis)),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
colnames(mtx.to.plot.val) <- plyr::mapvalues(gsub(" Sig ", " ", colnames(mtx.to.plot.val)),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)

ht.dis <- Heatmap(mtx.to.plot.dis,
        height = nrow(mtx.to.plot.dis) * unit(0.08, "in"),
        left_annotation = rowAnnotation(
                Cohort = sig.loading.mtx.df.list[["collapsed"]]$cohort,
                col = list(Cohort = cohort.pal),
                show_annotation_name = FALSE
        ),
        top_annotation = columnAnnotation(
                df = sig.loading.var.to.plot.dis |> select(!contains("mean")),
                # `Patient loading variance` = anno_boxplot(sig.loading.var.mtx, outline = FALSE,
                #                                           gp = gpar(fill = "white")),
                # `Patient loading means` = anno_boxplot(mtx.to.plot, outline = FALSE,
                #                                        gp = gpar(fill = "white"))
                show_annotation_name = FALSE
        ),
        column_split = sig.split,
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_names_rot = 75,
        name = "Signature loading",
        col = viridisLite::viridis(100, option = "C")
)

ht.val <- Heatmap(mtx.to.plot.val,
        height = nrow(mtx.to.plot.val) * unit(0.08, "in"),
        left_annotation = rowAnnotation(
                Cohort = sig.loading.mtx.df.list[["collapsed-scored-validation"]]$cohort,
                col = list(Cohort = cohort.pal),
                show_annotation_name = FALSE
        ),
        top_annotation = columnAnnotation(
                df = sig.loading.var.to.plot.val |> select(!contains("mean")),
                # `Patient loading variance` = anno_boxplot(sig.loading.var.mtx, outline = FALSE,
                #                                           gp = gpar(fill = "white")),
                # `Patient loading means` = anno_boxplot(mtx.to.plot, outline = FALSE,
                #                                        gp = gpar(fill = "white"))
                show_annotation_name = FALSE
        ),
        column_split = sig.split,
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_names_rot = 75,
        name = "Signature loading",
        col = viridisLite::viridis(100, option = "C")
)

ht <- Heatmap(rbind(mtx.to.plot.dis, mtx.to.plot.val),
        height = (nrow(mtx.to.plot.dis) + nrow(mtx.to.plot.val)) * unit(0.04, "in"),
        left_annotation = rowAnnotation(
                Cohort = c(sig.loading.mtx.df.list[["collapsed"]]$cohort, sig.loading.mtx.df.list[["collapsed-scored-validation"]]$cohort),
                col = list(Cohort = cohort.pal),
                show_annotation_name = FALSE,
                annotation_legend_param = list(direction = "horizontal")
        ),
        top_annotation = columnAnnotation(
                df = sig.loading.var.to.plot |> select(!contains("mean"), -signature),
                # `Patient loading variance` = anno_boxplot(sig.loading.var.mtx, outline = FALSE,
                #                                           gp = gpar(fill = "white")),
                # `Patient loading means` = anno_boxplot(mtx.to.plot, outline = FALSE,
                #                                        gp = gpar(fill = "white"))
                show_annotation_name = FALSE,
                annotation_legend_param = list(direction = "horizontal")
        ),
        heatmap_legend_param = list(direction = "horizontal"),
        column_split = sig.split,
        row_split = c(rep("Discovery", nrow(mtx.to.plot.dis)), rep("Validation", nrow(mtx.to.plot.val))),
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
        column_names_rot = 90,
        name = "Program/Gene\nloading",
        col = viridisLite::viridis(100, option = "C")
)

png(snakemake@output[['figure2_c']], width = snakemake@params[['figure2_c_width']], height = snakemake@params[['figure2_c_height']], units = "in", res = 360)
(ht.dis %v% ht.val) |>
        draw(
                merge_legends = TRUE
        )
dev.off()

print("Patient signature profile plot successfully created")

# panel D
gene.loading.top.gene.df.to.plot.list <- readRDS(snakemake@input[['selected_sig_top_gene_loading_mtx_list']])

gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.to.plot.list$`pancreatic epithelial cell`
ht.epi <- Heatmap(gene.loading.top.gene.df.to.plot,
                  name = "Gene loading",
                  #width = unit(snakemake@params[['figure2_d_width']] - 3, "in"),
                  column_labels = plyr::mapvalues(gsub(" Rep ", " ", colnames(gene.loading.top.gene.df.to.plot)),
                                               from = sig.interpt$signature,
                                               to = sig.interpt$`short interpretation`,
                                               warn_missing = FALSE),
                  #row_labels = c("Acinar cell", "General drug sensitivity", "Basal A / EMT", " Heatshock response", "Classical A/B", " Proliferation/SN38 sensitivity", "Ductal cell", "ZFAS1/P4HA1/EIF4A2", "Mitochondria metabolism", "Basal B - COL17A1", "Classical A", "Basal B - IL32"),
                  row_names_side = "left",
                  row_dend_side = "right",
                  col = colorRamp2(seq(from = 0, to = 0.7, length.out = 100), viridisLite::viridis(100, option = "C")),
                  show_heatmap_legend = FALSE)

gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.to.plot.list$`fibroblast`
ht.fibro <- Heatmap(gene.loading.top.gene.df.to.plot,
                    name = "Gene loading",
                    #width = unit(snakemake@params[['figure2_d_width']] - 3, "in"),
                    column_labels = plyr::mapvalues(gsub(" Rep ", " ", colnames(gene.loading.top.gene.df.to.plot)),
                                                 from = sig.interpt$signature,
                                                 to = sig.interpt$`short interpretation`,
                                                 warn_missing = FALSE),
                    #row_labels = c("Acinar cell", "General drug sensitivity", "Basal A / EMT", " Heatshock response", "Classical A/B", " Proliferation/SN38 sensitivity", "Ductal cell", "ZFAS1/P4HA1/EIF4A2", "Mitochondria metabolism", "Basal B - COL17A1", "Classical A", "Basal B - IL32"),
                    row_names_side = "left",
                    row_dend_side = "right",
                    col = colorRamp2(seq(from = 0, to = 0.7, length.out = 100), viridisLite::viridis(100, option = "C")),
                    show_heatmap_legend = FALSE)

gene.loading.top.gene.df.to.plot <- gene.loading.top.gene.df.to.plot.list$`CD8-positive, alpha-beta T cell`
ht.cd8 <- Heatmap(gene.loading.top.gene.df.to.plot,
                  name = "Gene loading",
                  #width = unit(snakemake@params[['figure2_d_width']] - 3, "in"),
                  column_labels = plyr::mapvalues(gsub(" Rep ", " ", colnames(gene.loading.top.gene.df.to.plot)),
                                               from = sig.interpt$signature,
                                               to = sig.interpt$`short interpretation`,
                                               warn_missing = FALSE),
                  #row_labels = c("Acinar cell", "General drug sensitivity", "Basal A / EMT", " Heatshock response", "Classical A/B", " Proliferation/SN38 sensitivity", "Ductal cell", "ZFAS1/P4HA1/EIF4A2", "Mitochondria metabolism", "Basal B - COL17A1", "Classical A", "Basal B - IL32"),
                  row_names_side = "left",
                  row_dend_side = "right",
                  col = colorRamp2(seq(from = 0, to = 0.7, length.out = 100), viridisLite::viridis(100, option = "C")),
                  show_heatmap_legend = FALSE)

png(snakemake@output[['figure2_d']], width = snakemake@params[['figure2_d_width']], height = snakemake@params[['figure2_d_height']], units = "in", res = 360)
ggarrange(as.grob(ht.epi), as.grob(ht.fibro), as.grob(ht.cd8),
          nrow = 1)
dev.off()

print("Selected signature top marker plot successfully created")

print("Gene loading correlation and distance plot successfully created")

# plot Figure 2
design <- "
AAAAABBB
CCCCCDDD
CCCCCDDD
CCCCCDDD
CCCCCDDD
"

ht.grob <- grid.grabExpr(draw(ht, merge_legend = TRUE))
ht.dist.corr.grob = grid.grabExpr(draw(ht.dist.corr, heatmap_legend_list = list(
        Legend(title = "Correlation", col_fun = col.corr),
        Legend(title = "Minkowski\nDistance", col_fun = col.dist)#,
        #Legend(title = "Correlation - Distance", col_fun = col.corrminusdist)
)))

ht.dist.corr.grob <- grid.grabExpr(draw(ht.coor + ht.dist, ht_gap = unit(-91.5, "mm")))

pdf(file = snakemake@output[["figure2_pdf"]], width = snakemake@params[["figure2_width"]], height = snakemake@params[["figure2_height"]])
p.signum + ht.dist.corr.grob + ht.grob + ggarrange(as.grob(ht.epi), as.grob(ht.fibro), as.grob(ht.cd8), nrow = 1) + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 18))
dev.off()

print("Figure 2 PDF successfully created")

png(file = snakemake@output[["figure2_png"]], width = snakemake@params[["figure2_width"]], height = snakemake@params[["figure2_height"]], units = "in", res = 360)
p.signum + ht.dist.corr.grob + ht.grob + ggarrange(as.grob(ht.epi), as.grob(ht.fibro), as.grob(ht.cd8), nrow = 1) + plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 18))
dev.off()

print("Figure 2 PNG successfully created")