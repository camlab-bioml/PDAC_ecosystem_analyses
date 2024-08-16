suppressPackageStartupMessages({
        library(tidyverse)
        library(ComplexHeatmap)
	library(ggplot2)
	library(ggpubr)
	library(ggsci)
        library(ggrepel)
	library(viridisLite)
	library(circlize)
        library(patchwork)
})

cell_type_rename <- read_csv(snakemake@input[["cell_type_rename"]])

sig.interpt <- readxl::read_xlsx(snakemake@input[["sig_interpretation"]])

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
        labs(
                color = "Same celltype",
                x = "Co-occurrence (Discovery)",
                y = "Co-occurrence (Validation)"
        ) +
        theme_pubr() +
        theme(
                legend.title = element_blank(),
                #legend.position = "bottom",
                legend.text = element_text(size = 14),
                axis.title.x = element_text(face = "bold"),
                axis.title.y = element_text(face = "bold"),
                text = element_text(face = "italic")
        ) +
        #guides(color = guide_legend(override.aes = list(size = 5, alpha = 1, shape = 16))) +
        guides(color = guide_legend(size = 12, nrow = 2, byrow = TRUE)) +
        stat_cor(aes(#label = after_stat(paste0("italic(R)~`=`~", r, "~`,`~italic(p)~`=`~", p))
                     #label = paste0(after_stat(r.label), "*','~", after_stat(p.label))
                     label = paste(after_stat(r.label), after_stat(p.label), sep = "*`,`~")
                     ))

ggsave(snakemake@output[["patient_profiles_correlation_comparison_overall"]], device = "png", width = 4.5, height = 5, units = "in", dpi = 360, bg = "white")

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

# plot signature co-occurrence in discovery vs. validation for each signature
dfc.list <- readRDS(snakemake@input[["signature_correlation_comparison_data_frame_list"]])

## make sure some signatures with same inpterpretation don't get stacked on each other
dfc.list <- lapply(dfc.list, function(dfc.sig) {
  dfc.sig$celltype <- plyr::mapvalues(dfc.sig$celltype, from = cell_type_rename$old_name, cell_type_rename$new_name)
  dfc.sig$signature <- plyr::mapvalues(dfc.sig$signature, 
                                         from = sig.interpt$signature, 
                                         to = sig.interpt$`short interpretation`,
                                       warn_missing = FALSE)
  dfc.sig |> 
    filter(!grepl("Ambient RNA|^MALAT1/NEAT1$", signature)) |>
    group_by(signature) |>
    mutate(unique_signature = paste0(signature, "_-_", row_number())) |>
    ungroup() |>
    arrange(desc(correlation_discovery))
})

## draw some bar/scatter plots
individual_plot_dir <- snakemake@params[["sig_cooccurring_plot_dir"]]

lapply(names(dfc.list), function(sig) {
  if(nrow(dfc.list[[sig]]) == 0) {
    return()
  }
  dfc.sig <- dfc.list[[sig]] |>
    mutate(sig.of.interest = ifelse(((abs(correlation_discovery) > 0.3 | abs(correlation_validation) > 0.3) & 
                                       abs(correlation_discovery - correlation_validation) <= 0.2), signature, NA))
  plot.title <- sig.interpt |> filter(signature == sig) |> pull(`short interpretation`)
#   plot.title <- paste0(plyr::mapvalues(gsub(" [0-9]$| [0-9][0-9]$", "", sig), 
#                                        from = cell_type_rename$old_name, 
#                                        to = cell_type_rename$new_name, warn_missing = FALSE), ": ", plot.title)
  plot.title <- paste0("Co-occurrence with ", plot.title)
  
  p1 <- ggplot(dfc.sig, aes(x = reorder(unique_signature, -correlation_discovery), y = correlation_discovery)) +
    geom_bar(stat = "identity", aes(fill = celltype, alpha = `validation confidence`)) +
    scale_fill_manual(values = celltype_pal_to_use) +
    labs(title = plot.title, x = NULL, y = "Co-occurrence in discovery") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(axis.text.x = element_blank(),
          plot.margin = margin(l = 0 + 200))
  
  p2 <- ggplot(dfc.sig, aes(x = reorder(unique_signature, -correlation_discovery), y = correlation_validation)) +
    geom_bar(stat = "identity", aes(fill = celltype, alpha = `validation confidence`)) +
    scale_x_discrete(labels = dfc.sig$signature) + 
    scale_fill_manual(values = celltype_pal_to_use) +
    labs(x = NULL, y = "Co-occurrence in validation") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(l = 0 + 200))
  p1 / p2 + plot_layout(guides = "collect") &
    theme(legend.position='top')
  ggsave(paste0(individual_plot_dir, sig, "-cooccur-bar.png"), width = 18, height = 15, units = "in", dpi = 360)
  
  ggplot(dfc.sig, aes(x = reorder(unique_signature, -correlation_discovery), y = correlation_discovery)) +
    geom_bar(stat = "identity", aes(alpha = `validation confidence`)) +
    scale_x_discrete(labels = dfc.sig$signature) + 
    labs(title = plot.title, x = NULL, y = "Co-occurrence in discovery") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(l = 0 + 200))
  ggsave(paste0(individual_plot_dir, sig,"-cooccur-bar-no-color.png"), width = 18, height = 10, units = "in", dpi = 360)
  
  ggplot(dfc.sig, aes(x = correlation_discovery, y = correlation_validation, label = sig.of.interest, fill = celltype, color = celltype)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_point(aes(fill = celltype, color = celltype), size = 4) +
    geom_label_repel(show.legend = F, max.overlaps = 15, alpha = 0.8, colour = "black", 
                     arrow = arrow(type = "closed", angle = 20, length = unit(0.08, "inches"))) +
    scale_fill_manual(values = celltype_pal_to_use) + 
    scale_color_manual(values = celltype_pal_to_use) +
    labs(title = plot.title,
         x = "Correlation in Discovery",
         y = "Correlation in Validation",
         color = "Cell type") + 
    theme_pubr() +
    theme(axis.title = element_text(size = 16, face = "bold"),
          title = element_text(face = "bold")) + 
    guides(fill = "none")
  ggsave(paste0(individual_plot_dir, sig, "-cooccur-scatter.png"), width = 7, height = 6, units = "in", dpi = 360)
})

# plot general signature co-occurrence agreement between discovery an validation
dis.val.agree <- read_tsv(snakemake@input[["signature_cooccurrence_agreement_data_frame"]])

dis.val.agree <- dis.val.agree |>
  mutate(celltype = plyr::mapvalues(celltype, from = cell_type_rename$old_name, to = cell_type_rename$new_name)) |>
  mutate(signature = plyr::mapvalues(signature, from = sig.interpt$signature, to = sig.interpt$`short interpretation`))

dis.val.agree <- dis.val.agree |> 
  filter(!grepl("Ambient RNA|^MALAT1/NEAT1$", signature)) |>
  group_by(signature) |>
  mutate(unique_signature = paste0(signature, "_-_", row_number())) |>
  ungroup() |>
  arrange(desc(dis_val_corr))

margin_spacer <- function(x) {
  # where x is the column in your dataset
  #left_length <- nchar(levels(factor(x)))[1]
  left_length <- nchar(x)[1]
  if (left_length > 8) {
    return((left_length - 8) * 4)
  }
  else
    return(0)
}

ggplot(dis.val.agree, aes(x = reorder(unique_signature, -dis_val_corr), y = dis_val_corr)) +
  geom_bar(stat = "identity", aes(fill = celltype)) +
  ylim(0, NA) +
  scale_x_discrete(labels = dis.val.agree$signature) + 
  scale_fill_manual(values = celltype_pal_to_use) +
  labs(x = NULL, y = "Correlation of co-occurrence between discovery and validation") +
  geom_text(aes(label = ifelse(p_value < 0.01, "*", "")), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        #plot.margin = margin(l = 0 + margin_spacer(dis.val.agree$signature))
        plot.margin = margin(l = 0 + 32))
ggsave(snakemake@output[["signature_cooccurrence_agreement_plot"]], width = 18, height =8, units = "in", dpi = 360)