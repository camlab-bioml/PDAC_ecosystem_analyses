suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(scales)
  library(BiocParallel)
  library(sjstats)
  library(dittoSeq)
  library(ggplotify)
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
        "Intra-patient heterogeneity \n Discovery",
        "Inter-patient heterogeneity \n Discovery (mean)",
        "Inter-patient heterogeneity \n Discovery")
names(sig.loading.var.to.plot.val) <- c(
        "Intra-patient heterogeneity \n Validation",
        "Inter-patient heterogeneity \n Validation (mean)",
        "Inter-patient heterogeneity \n Validation")

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
        select(-ambient.sigs) |>
        reframe(across(where(is.numeric), rescale)) |>
        as.matrix()
mtx.to.plot.val <- sig.loading.mtx.df.list[["collapsed-scored-validation"]] |>
        select(-ambient.sigs) |>
        reframe(across(where(is.numeric), rescale)) |>
        as.matrix()

# colnames(mtx.to.plot.dis) <- paste0(sig.split, " Sig ", str_split(colnames(mtx.to.plot.dis), " Sig ", simplify = TRUE)[,2])
# colnames(mtx.to.plot.val) <- paste0(sig.split, " Sig ", str_split(colnames(mtx.to.plot.val), " Sig ", simplify = TRUE)[,2])

colnames(mtx.to.plot.dis) <- plyr::mapvalues(gsub(" Sig ", " ", colnames(mtx.to.plot.dis)),
        from = sig.interpt$signature,
        to = sig.interpt$`intepretation (change interpretation to have main marker for the unclear sigs, and put the interpretation a sub interpretation)`
)
colnames(mtx.to.plot.val) <- plyr::mapvalues(gsub(" Sig ", " ", colnames(mtx.to.plot.val)),
        from = sig.interpt$signature,
        to = sig.interpt$`intepretation (change interpretation to have main marker for the unclear sigs, and put the interpretation a sub interpretation)`
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

png(snakemake@output[['figure2_b']], width = snakemake@params[['figure2_b_width']], height = snakemake@params[['figure2_b_height']], units = "in", res = 321)
(ht.dis %v% ht.val) |>
        draw(
                merge_legends = TRUE
        )
dev.off()