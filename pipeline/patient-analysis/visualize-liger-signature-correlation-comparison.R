suppressPackageStartupMessages({
        library(tidyverse)
        library(ComplexHeatmap)
	library(ggplot2)
	library(ggpubr)
	library(ggsci)
	library(viridisLite)
	library(circlize)
})

cell_type_rename <- read_csv(snakemake@input[["cell_type_rename"]])

# signatures to remove
ambient.sigs <- read_csv(snakemake@input[["ambient_sigs"]])$signature

ambient.sigs <- gsub(" Sig ", " ", ambient.sigs)

# load cohort color palette
cohort_pal <- readRDS(snakemake@input[["cohort_pal"]])

# load celltype color palette
celltype_pal <- readRDS(snakemake@input[["celltype_pal"]])

celltype_pal$Cell_type <- str_split(celltype_pal$Cell_type_dis, " \\(", simplify = TRUE)[, 1]

celltype_pal_to_use <- celltype_pal$color_dis
names(celltype_pal_to_use) <- celltype_pal$Cell_type

# Load the data
dfc <- read_tsv(snakemake@input[["patient_profiles_correlation_data_frame"]])

## Discovery cohort
patient.profiles.df3 <- read_tsv(snakemake@input[["patient_profiles_collapsed_cleaned_scaled"]])
patient.profiles.df3 <- patient.profiles.df3 |> select(-cohort)

patient.profiles.mat <- select(patient.profiles.df3, -sample) |>
        as.matrix()
rownames(patient.profiles.mat) <- patient.profiles.df3$sample

patient.profiles.meta <- read_tsv(snakemake@input[["patient_profiles_collapsed_meta"]])
patient.profiles.mat <- patient.profiles.mat[patient.profiles.meta$sample, ]

## Validation cohort
patient.validation.profiles.df3 <- read_tsv(snakemake@input[["patient_profiles_collapsed_scored_val_cleaned_scaled"]])
patient.validation.profiles.df3 <- patient.validation.profiles.df3 |> select(-cohort)

patient.validation.profiles.mat <- select(patient.validation.profiles.df3, -sample) |>
        as.matrix()
rownames(patient.validation.profiles.mat) <- patient.validation.profiles.df3$sample

patient.validation.profiles.meta <- read_tsv(snakemake@input[["patient_profiles_collapsed_scored_val_meta"]])
patient.validation.profiles.mat <- patient.validation.profiles.mat[patient.validation.profiles.meta$sample, ]


# heatmap
## Discovery cohort
set.seed(42L)
top_annot <- columnAnnotation(
  cohort = patient.profiles.meta$cohort,
  col = list(cohort = cohort_pal)
)
right_annot <- rowAnnotation(
  celltype = plyr::mapvalues(gsub(" Rep [0-9]$| Rep [0-9][0-9]$", "", colnames(patient.profiles.mat)),
			     from = cell_type_rename$old_name,
			     to = cell_type_rename$new_name),
  col = list(celltype = celltype_pal_to_use)
)

patient.profiles.mat.for.plot <- patient.profiles.mat
patient.profiles.mat.for.plot[is.na(patient.profiles.mat.for.plot)] <- 0

png(filename = snakemake@output[["patient_profiles_collapsed_cleaned_sclaed_heatmap"]], width = 15, height = 10, units = "in", res = 321)
Heatmap(t(patient.profiles.mat.for.plot),
	col = colorRamp2(seq(from = 0, to = 1, length.out = 20), viridis(20)),
        top_annotation = top_annot,
        right_annotation = right_annot)
dev.off()

## Validation cohort
set.seed(42L)
top_annot <- columnAnnotation(
  cohort = patient.validation.profiles.meta$cohort,
  col = list(cohort = cohort_pal)
)
right_annot <- rowAnnotation(
  celltype = plyr::mapvalues(gsub(" RepVal [0-9]$| RepVal [0-9][0-9]$", "", colnames(patient.validation.profiles.mat)),
			     from = cell_type_rename$old_name,
			     to = cell_type_rename$new_name),
  col = list(celltype = celltype_pal_to_use)
)

patient.validation.profiles.mat.for.plot <- patient.validation.profiles.mat
patient.validation.profiles.mat.for.plot[is.na(patient.validation.profiles.mat.for.plot)] <- 0

png(filename = snakemake@output[["patient_profiles_collapsed_scored_val_cleaned_scaled_heatmap"]], width = 15, height = 10, units = "in", res = 321)
Heatmap(t(patient.validation.profiles.mat.for.plot),
	col = colorRamp2(seq(from = 0, to = 1, length.out = 20), viridis(20)),
	top_annotation = top_annot,
	right_annotation = right_annot)
dev.off()

# see how these match
## heatmap
colnames(patient.profiles.mat) <- gsub(" Rep", "", colnames(patient.profiles.mat))
colnames(patient.validation.profiles.mat) <- gsub(" RepVal", "", colnames(patient.validation.profiles.mat))
common.sigs <- intersect(colnames(patient.profiles.mat), colnames(patient.validation.profiles.mat))

png(filename = snakemake@output[["patient_profiles_combined_cleaned_scaled_heatmap"]], width = 20, height = 10, units = "in", res = 321)
Heatmap(t(patient.profiles.mat)[common.sigs, ], name = "validated/discovery",
	col = colorRamp2(seq(from = 0, to = 1, length.out = 20), viridis(20))) +
        Heatmap(t(patient.validation.profiles.mat)[common.sigs, ], name = "validation",
		col = colorRamp2(seq(from = 0, to = 1, length.out = 20), viridis(20)))
dev.off()

# Are signature correlations consistent between cohorts?
## remove ambient signatures
dfc.for.plot <- dfc
dfc.for.plot <- dfc.for.plot |>
        filter(!(signature_1 %in% ambient.sigs)) |>
        filter(!(signature_2 %in% ambient.sigs))

## update cell type labels for plotting
dfc.for.plot$cell_type_1 <- plyr::mapvalues(dfc.for.plot$cell_type_1,
                                            from = cell_type_rename$old_name,
                                            to = cell_type_rename$new_name)
dfc.for.plot$cell_type_2 <- plyr::mapvalues(dfc.for.plot$cell_type_2,
                                            from = cell_type_rename$old_name,
                                            to = cell_type_rename$new_name)

## tidy up for plotting
dfc.for.plot <- dfc.for.plot |>
        mutate(same_cell_type_for_plot = ifelse(same_cell_type,
                "Co-occurrence within cell type",
                "Co-occurrence between cell types"
        ))

## plot overall correlations
ggplot(dfc.for.plot, aes(x = correlation_discovery, y = correlation_validation)) +
        geom_point(aes(color = same_cell_type_for_plot), alpha = 0.5) +
        # facet_wrap(~ cell_type_1, scales = "free", nrow = 2) +
        geom_smooth(method = "lm", colour = "grey30") +
        stat_cor() +
        labs(
                color = "Same celltype",
                x = "Signature co-occurrence in discovery",
                y = "Signature co-occurrence in validation"
        ) +
        theme_pubr() +
        theme(
                legend.title = element_blank(),
                axis.title.x = element_text(face = "bold"),
                axis.title.y = element_text(face = "bold")
        ) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave(snakemake@output[["patient_profiles_correlation_comparison_overall"]], device = "png", width = 5, height = 5.5, units = "in", dpi = 321, bg = "white")

# plot intra-cell type correlations
filter(dfc.for.plot, same_cell_type) |>
        ggplot(aes(x = correlation_discovery, y = correlation_validation)) +
        geom_point(aes(colour = cell_type_1)) +
        scale_color_manual(values = celltype_pal_to_use) +
        facet_wrap(~cell_type_1, scales = "free", nrow = 2) +
        geom_smooth(method = "lm") +
        stat_cor() +
        labs(
                colour = "Celltype",
                x = "Signature co-occurrence in discovery",
                y = "Signature co-occurrence in validation"
        ) +
        theme_pubr() +
        theme(axis.title = element_text(face = "bold"))

ggsave(snakemake@output[["patient_profiles_correlation_comparison_intra_celltype"]], device = "png", width = 10, height = 5, units = "in", dpi = 321, bg = "white")

# plot inter-cell type correlations
dfc.for.plot_list <- list()

for (ct in union(unique(dfc.for.plot$cell_type_1), unique(dfc.for.plot$cell_type_2))) {
        dfc.for.plot_list[[ct]] <- dfc.for.plot %>%
                filter(!same_cell_type) %>%
                filter(cell_type_1 == ct | cell_type_2 == ct) |>
                mutate(facet_cell_type = ct) |>
                mutate(cell_type = ifelse(cell_type_1 == ct, cell_type_2, cell_type_1))
}
rm(ct)

Reduce(rbind, dfc.for.plot_list) |>
        ggplot(aes(x = correlation_discovery, y = correlation_validation)) +
        geom_point(aes(colour = cell_type)) +
        scale_color_manual(values = celltype_pal_to_use) +
        facet_wrap(~facet_cell_type, scales = "free", nrow = 2) +
        geom_smooth(method = "lm") +
        stat_cor() +
        labs(
                colour = "Celltype",
                x = "Signature co-occurrence in discovery",
                y = "Signature co-occurrence in validation"
        ) +
        theme_pubr() +
        theme(axis.title = element_text(face = "bold"))

ggsave(snakemake@output[["patient_profiles_correlation_comparison_inter_celltype"]], device = "png", width = 10, height = 5, units = "in", dpi = 321, bg = "white")