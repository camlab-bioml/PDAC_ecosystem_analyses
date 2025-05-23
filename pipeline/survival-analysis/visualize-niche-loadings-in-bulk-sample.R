suppressPackageStartupMessages({
        library(magrittr)
        library(tidyverse)
        library(scales)
        library(sjstats)
        library(ggplotify)
        library(ggsci)
        library(ggpubr)
        library(ggpmisc)
        library(cowplot)
        library(circlize)
        library(ComplexHeatmap)
        library(patchwork)
        library(magick)
        library(see)
})

# load data
chan.niche.loadings.gsva <- readRDS(snakemake@input[["gsva_es_pancurx_niche"]])
chan.niche.loadings.exprs <- readRDS(snakemake@input[["sig_exprs_mean_pancurx_niche"]])
chan.clin.data <- read_tsv(snakemake@input[["tidy_pancurx_clin_data"]])

tcga.niche.loadings.gsva <- readRDS(snakemake@input[["gsva_es_paad_niche"]])
tcga.niche.loadings.exprs <- readRDS(snakemake@input[["sig_exprs_mean_paad_niche"]])
tcga.clin.data <- read_tsv(snakemake@input[["tidy_paad_clin_data"]])

print("niche loading data frames")
print(head(chan.niche.loadings.exprs))
print(nrow(chan.niche.loadings.exprs))
print(head(tcga.niche.loadings.exprs))
print(nrow(tcga.niche.loadings.exprs))

# tidy up
## TCGA
tcga.niche.loadings.exprs <- tcga.niche.loadings.exprs |>
    column_to_rownames("bcr_patient_barcode") |>
    select(-Niche_1) |>
    as.matrix()
colnames(tcga.niche.loadings.exprs) <- c("Classical", "Basal-like", "Immune act.")

tcga.niche.loadings.gsva <- tcga.niche.loadings.gsva |>
    column_to_rownames("bcr_patient_barcode") |>
    select(-Niche_1) |>
    as.matrix()
colnames(tcga.niche.loadings.gsva) <- c("Classical", "Basal-like", "Immune act.")

tcga.sample.id.tmp <- tcga.niche.loadings.gsva |> as.data.frame() |>
    rownames_to_column("bcr_patient_barcode") |>
    select(bcr_patient_barcode) |>
    mutate(bcr_patient_barcode = str_sub(bcr_patient_barcode, 1, 12))

print(head(tcga.sample.id.tmp))

tcga.clin.data <- left_join(tcga.sample.id.tmp, tcga.clin.data, by = "bcr_patient_barcode") |>
    column_to_rownames("bcr_patient_barcode") |>
    select(ajcc_pathologic_tumor_stage, new_tumor_event_type, treatment_outcome_first_course)

tcga.clin.data <- tcga.clin.data |>
    mutate(ClinicalStaging = case_when(
	ajcc_pathologic_tumor_stage %in% c("Stage IA", "Stage IB", "Stage IIA", "Stage IIB") ~ "Early",
	ajcc_pathologic_tumor_stage %in% c("Stage III/IV") ~ "Late",
	TRUE ~ NA_character_
    )) |>
    mutate(NewTumorEvent = case_when(
	new_tumor_event_type %in% c("Locoregional Recurrence|Distant Metastasis") ~ "LR|DM",
	new_tumor_event_type %in% c("Distant Metastasis") ~ "DM",
	new_tumor_event_type %in% c("Locoregional Recurrence") ~ "LR",
	new_tumor_event_type %in% c("New Primary Tumor") ~ "NPT",
	TRUE ~ NA_character_
    )) |>
    mutate(TreatmentOutcome = case_when(
	treatment_outcome_first_course %in% c("Complete Remission/Response") ~ "CR",
	treatment_outcome_first_course %in% c("Partial Remission/Response") ~ "PR",
	treatment_outcome_first_course %in% c("Stable Disease") ~ "SD",
	treatment_outcome_first_course %in% c("Progressive Disease") ~ "PD",
	TRUE ~ NA_character_
    )) |>
    select(ClinicalStaging, NewTumorEvent, TreatmentOutcome)

print(head(tcga.clin.data))

row_ha_tcga <- rowAnnotation(
    df = tcga.clin.data
)

## Chan
chan.niche.loadings.exprs <- chan.niche.loadings.exprs |>
    column_to_rownames("sample_id") |>
    select(-Niche_1) |>
    as.matrix()
colnames(chan.niche.loadings.exprs) <- c("Classical", "Basal-like", "Immune act.")

chan.niche.loadings.gsva <- chan.niche.loadings.gsva |>
    column_to_rownames("sample_id") |>
    select(-Niche_1) |>
    as.matrix()
colnames(chan.niche.loadings.gsva) <- c("Classical", "Basal-like", "Immune act.")

chan.sample.id.tmp <- chan.niche.loadings.gsva |> as.data.frame() |>
    rownames_to_column("sample_id") |>
    select(sample_id) 

print(head(chan.sample.id.tmp))

chan.clin.data <- left_join(chan.sample.id.tmp, chan.clin.data, by = "sample_id") |>
    column_to_rownames("sample_id") |>
    select(ploidy, purity, ClinicalStaging, AdjuvantBest_Therapy_Response)

chan.clin.data <- chan.clin.data |>
    mutate(ClinicalStaging = case_when(
	ClinicalStaging %in% c("IA", "IB", "IB or IIB", "II", "IIA", "IIB", "IIA or IIB") ~ "Early",
	ClinicalStaging %in% c("III", "I-III") ~ "Late",
	ClinicalStaging %in% c("Unknown") ~ NA_character_
    )) |>
    mutate(AdjBestResponse = case_when(
	AdjuvantBest_Therapy_Response %in% c("NED") ~ "NED",
	AdjuvantBest_Therapy_Response %in% c("Disease Progression") ~ "PD",
	TRUE ~ NA_character_
    )) |>
    select(ClinicalStaging, AdjBestResponse, ploidy, purity)

print(head(chan.clin.data))

row_ha_chan <- rowAnnotation(
    df = chan.clin.data
)

# draw heatmap
ht.chan <- Heatmap(
	scale(chan.niche.loadings.gsva), 
	name = "Ecosystem loadings", 
	#col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
	show_row_names = FALSE, 
	left_annotation = row_ha_chan,
	show_column_names = TRUE, 
	column_names_rot = 90,
	column_title = "Ecosystem", 
	row_title = "Patient", 
	cluster_rows = TRUE, 
	cluster_columns = FALSE, 
	column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
	row_title_gp = gpar(fontsize = 12, fontface = "bold"), 
	heatmap_legend_param = list(
		title = "Ecosystem\nloading\nChan", 
		#at = seq(-1, 1, 0.5), 
		#labels = c("-1", "-0.5", "0", "0.5", "1"), 
		direction = "vertical", 
		#legend_width = unit(2, "cm"), 
		#legend_height = unit(0.5, "cm"), 
		legend_gp = gpar(fontsize = 10)
	),
	row_km = 4,
	border = TRUE
)

ht.tcga <- Heatmap(
	scale(tcga.niche.loadings.gsva), 
	name = "Ecosystem loadings", 
	#col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
	show_row_names = FALSE, 
	left_annotation = row_ha_tcga,
	show_column_names = TRUE, 
	column_names_rot = 90,
	column_title = "Ecosystem", 
	row_title = "Patient", 
	cluster_rows = TRUE, 
	cluster_columns = FALSE, 
	column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
	row_title_gp = gpar(fontsize = 12, fontface = "bold"), 
	heatmap_legend_param = list(
		title = "Ecosystem\nloading\nTCGA", 
		#at = seq(-1, 1, 0.5), 
		#labels = c("-1", "-0.5", "0", "0.5", "1"), 
		direction = "vertical", 
		#legend_width = unit(2, "cm"), 
		#legend_height = unit(0.5, "cm"), 
		legend_gp = gpar(fontsize = 10)
	),
	row_km = 4,
	border = TRUE
)

# draw scatter plot
print(colnames(as.data.frame(chan.niche.loadings.gsva)))
print(colnames(as.data.frame(tcga.niche.loadings.gsva)))

colnames(chan.niche.loadings.gsva) <- c("Classical", "Basal", "Immune")
colnames(tcga.niche.loadings.gsva) <- c("Classical", "Basal", "Immune")

sc.chan <- ggscatter(
        data = as.data.frame(chan.niche.loadings.gsva),
        x = "Classical",
        y = "Basal",
        color = "Immune",
        shape = 19,
        size = 2.5,
        alpha = 0.6,
        xlab = "Classical",
        ylab = "Basal-like",
        title = "GSVA loadings - Chan",
        legend.title = "Immune act.",
        legend.position = "right"
) + 
scale_color_gradient(low = "blue", high = "red")

sc.tcga <- ggscatter(
        data = as.data.frame(tcga.niche.loadings.gsva),
        x = "Classical",
        y = "Basal",
        color = "Immune",
        shape = 19,
        size = 2.5,
        alpha = 0.6,
        xlab = "Classical",
        ylab = "Basal-like",
        title = "GSVA loadings - TCGA",
        legend.title = "Immune act.",
        legend.position = "right"
) + 
scale_color_gradient(low = "blue", high = "red")

# save plots
png(snakemake@output[["figure_niche_loadings_paad"]], width = 4, height = 6, units = "in", res = 600)
draw(ht.tcga)
dev.off()

png(snakemake@output[["figure_niche_loadings_pancurx"]], width = 4, height = 6, units = "in", res = 600)
draw(ht.chan)
dev.off()

ggsave(snakemake@output[["figure_scatter_niche_loadings_pancurx"]], sc.chan, width = 6, height = 6, dpi = 600)
ggsave(snakemake@output[["figure_scatter_niche_loadings_paad"]], sc.tcga, width = 6, height = 6, dpi = 600)