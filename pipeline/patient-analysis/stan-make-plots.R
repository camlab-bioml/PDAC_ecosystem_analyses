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
  library(ComplexHeatmap)
  library(patchwork)
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

# sample numeric encoding used in stan
samples.encodings <- read_tsv(snakemake@input[["samples_encodings"]])

# load cohort color palette
cohort_pal <- readRDS(snakemake@input[["cohort_pal"]])

# load celltype color palette
celltype_pal <- readRDS(snakemake@input[["celltype_pal"]])

celltype_pal$Cell_type <- str_split(celltype_pal$Cell_type_dis, " \\(", simplify = TRUE)[, 1]

celltype_pal_to_use <- celltype_pal$color_dis
names(celltype_pal_to_use) <- celltype_pal$Cell_type

# load stan data and other parameters
celltypes <- snakemake@params[["celltypes"]]
number.of.niches <- snakemake@params[["number_of_niches"]]
nIter <- snakemake@params[["nIter"]]

niche.sig.loading.init.value <- readRDS(snakemake@input[["microenvironment_niche_factors_init_value"]])
niche.loading.init.value <- readRDS(snakemake@input[["niche_factor_loadings_init_value"]])

infered.sig.loading <- readRDS(snakemake@input[['patient_specific_modelled_mu']])
niche.sig.loading <- readRDS(snakemake@input[['microenvironment_niche_factors']])
niche.loading <- readRDS(snakemake@input[['niche_factor_loadings']])
cov.i.list <- readRDS(snakemake@input[["intrinsic_covariance_matrices"]])

# plot patient specific modelled mu - infered signature loading from stan
column.split <- plyr::mapvalues(gsub(" [0-9]$| [0-9][0-9]$", "", colnames(infered.sig.loading)),
                                from = cell_type_rename$old_name,
                                to = cell_type_rename$new_name)

colnames(infered.sig.loading) <- plyr::mapvalues(colnames(infered.sig.loading),
  from = sig.interpt$signature,
  to = sig.interpt$`short interpretation`
)

png(snakemake@output[["patient_specific_modelled_mu_plot"]], width = 17, height = 16, units = "in", res = 321)
Heatmap(infered.sig.loading,
        height = nrow(infered.sig.loading) * unit(0.07, "in"),
        name = "Estimated\nsignature\nloading",
        top_annotation = columnAnnotation(
          celltype = column.split,
          col = list(celltype = celltype_pal_to_use),
          show_annotation_name = FALSE
        ),
        column_names_rot = 90,
        col = viridisLite::viridis(100, option = "C")
        ) |>
  draw(column_title = paste0("Patient specific modelled mu", " - ", snakemake@wildcards[["condition"]], " - ", number.of.niches, " niches, ", nIter, " iterations"))
dev.off()

# plot microenvironment niche factors - infered signature loading in niches from stan
column.split <- plyr::mapvalues(gsub(" [0-9]$| [0-9][0-9]$", "", colnames(niche.sig.loading)),
  from = cell_type_rename$old_name,
  to = cell_type_rename$new_name
)

colnames(niche.sig.loading) <- plyr::mapvalues(colnames(niche.sig.loading),
  from = sig.interpt$signature,
  to = sig.interpt$`short interpretation`
)

png(snakemake@output[["microenvironment_niche_factors_plot"]], width = 17, height = 13, units = "in", res = 321)
Heatmap(niche.sig.loading,
        height = nrow(niche.sig.loading) * unit(0.4, "in"),
        name = "Signature\nloading",
        top_annotation = columnAnnotation(
          celltype = column.split,
          col = list(celltype = celltype_pal_to_use),
          show_annotation_name = FALSE
        ),
        column_names_rot = 90,
        col = viridisLite::viridis(100, option = "C")) |> 
  draw(column_title = paste0("Microenvironment Niche Factors", " - ", snakemake@wildcards[["condition"]], " - ", number.of.niches, " niches, ", nIter, " iterations"))
dev.off()

## plot microenvironment niche factors - initialized values
colnames(niche.sig.loading.init.value) <- plyr::mapvalues(colnames(niche.sig.loading.init.value),
  from = sig.interpt$signature,
  to = sig.interpt$`short interpretation`
)

png(snakemake@output[["microenvironment_niche_factors_init_value_plot"]], width = 17, height = 13, units = "in", res = 321)
Heatmap(niche.sig.loading.init.value,
        height = nrow(niche.sig.loading.init.value) * unit(0.4, "in"),
        name = "Signature\nloading",
        top_annotation = columnAnnotation(
          celltype = column.split,
          col = list(celltype = celltype_pal_to_use),
          show_annotation_name = FALSE
        ),
        column_names_rot = 90,
        col = viridisLite::viridis(100, option = "C")) |> 
  draw(column_title = paste0("Initialized Microenvironment Niche Factors", " - ", snakemake@wildcards[["condition"]], " - ", number.of.niches, " niches, ", nIter, " iterations"))
dev.off()

# plot niche factor loadings - infered niche loading in samples from stan
png(snakemake@output[["niche_factor_loadings_plot"]], width = 7, height = 10, units = "in", res = 321)
Heatmap(niche.loading,
        name = "Niche\nloading",
        right_annotation = rowAnnotation(
          cohort = plyr::mapvalues(rownames(niche.loading),
                                   from = samples.encodings$sample,
                                   to = samples.encodings$cohort,
                                   warn_missing = FALSE),
          col = list(cohort = cohort_pal),
          show_annotation_name = FALSE
        ),
        show_row_names = FALSE,
        col = viridisLite::viridis(100, option = "C")) |> 
  draw(column_title = paste0("Niche loadings", " - ", snakemake@wildcards[["condition"]], " - ", number.of.niches, " niches, ", nIter, " iterations"))
dev.off()

# plot niche factor loadings - initialized values
png(snakemake@output[["niche_factor_loadings_init_value_plot"]], width = 7, height = 10, units = "in", res = 321)
Heatmap(niche.loading.init.value,
        name = "Niche\nloading",
        right_annotation = rowAnnotation(
          cohort = plyr::mapvalues(rownames(niche.loading.init.value),
                                   from = samples.encodings$sample,
                                   to = samples.encodings$cohort,
                                   warn_missing = FALSE),
          col = list(cohort = cohort_pal),
          show_annotation_name = FALSE
        ),
        show_row_names = FALSE,
        col = viridisLite::viridis(100, option = "C")) |> 
  draw(column_title = paste0("Initialized Niche loadings", " - ", snakemake@wildcards[["condition"]], " - ", number.of.niches, " niches, ", nIter, " iterations"))
dev.off()

# plot intrinsic covariance matrices - infered intrinsic covariance in samples from stan
cov.i.ht.list <- lapply(names(cov.i.list), function(ct) {
  cov.i <- cov.i.list[[ct]]
  cov.i.ht <- Heatmap(cov.i,
    name = paste0(snakemake@wildcards[["condition"]], "\nintrinsic covariance"),
    show_row_names = FALSE,
    # show_column_names = FALSE,
    # cluster_rows = FALSE, cluster_columns = FALSE,
    # column_title = "cov(i)", column_title_gp = gpar(fontsize = 20),
    # row_title = "cov(i)", row_title_gp = gpar(fontsize = 20),
    # heatmap_legend_param = list(title = "cov(i)", title_gp = gpar(fontsize = 20)),
    col = circlize::colorRamp2(c(min(cov.i), 0, max(cov.i)), c("blue", "white", "red"))
  )
  cov.i.ht
})

cov.i.ht.list <- lapply(cov.i.ht.list, function(cov.i.ht) draw(cov.i.ht) |> grid.grabExpr())

png(snakemake@output[["intrinsic_covariance_matrices_plot"]], width = 20, height = 10, units = "in", res = 321)
# ComplexHeatmap::draw(ht.list, merge_legends = TRUE)
# plot_grid(cov.i.ht.list, nrow = 2)
wrap_plots(cov.i.ht.list, nrow = 2)
dev.off()