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
  library(survival)
  library(survminer)
  library(circlize)
  library(patchwork)
  library(magick)
  library(see)
  library(ggrepel)
  #library(forestplot)
})

celltype_rename <- read_csv(snakemake@input[['cell_type_rename']])

# load signature interpretation
sig.interpt <- readxl::read_xlsx(snakemake@input[['sig_interpretation']])

# signatures to remove
ambient.sigs <- read_csv(snakemake@input[['ambient_sigs']])$signature

# load celltype color palette
celltype_pal <- readRDS(snakemake@input[["celltype_pal"]])
celltype_pal$Cell_type <- str_split(celltype_pal$Cell_type_dis, " \\(", simplify = TRUE)[, 1]
celltype_pal_to_use <- celltype_pal$color_dis
names(celltype_pal_to_use) <- celltype_pal$Cell_type

# Panel E - supp
## tcga
clin.assoc.tcga <- readRDS(snakemake@input[["clin_assoc_paad"]]) |>
        mutate(celltype = str_split(signature, "\\.", simplify = T)[, 1]) |>
        mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
clin.assoc.tcga$celltype <- plyr::mapvalues(clin.assoc.tcga$celltype,
        from = celltype_rename$old_name,
        to = celltype_rename$new_name
)
clin.assoc.tcga$signature <- plyr::mapvalues(gsub(" Rep ", " ", clin.assoc.tcga$signature),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
clin.assoc.tcga <- clin.assoc.tcga |>
        filter(!grepl("Ambient RNA|MALAT1|HSP|Proliferation", signature)) |>
        column_to_rownames("signature") |>
        mutate(log_HR_tcga = `log(HR)`)

clin.assoc.tcga.pval <- readRDS(snakemake@input[["clin_assoc_paad_pval"]]) |> 
        mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
clin.assoc.tcga.pval$signature <- plyr::mapvalues(gsub(" Rep ", " ", clin.assoc.tcga.pval$signature),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
clin.assoc.tcga.pval <- clin.assoc.tcga.pval |>
        filter(!grepl("Ambient RNA|MALAT1|HSP|Proliferation", signature)) |>
        column_to_rownames("signature") |>
        mutate(log_HR_tcga = `log(HR)`)

## pancurx
clin.assoc.pancurx <- readRDS(snakemake@input[["clin_assoc_pancurx"]]) |>
        mutate(celltype = str_split(signature, "\\.", simplify = T)[, 1]) |>
        mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
clin.assoc.pancurx$celltype <- plyr::mapvalues(clin.assoc.pancurx$celltype,
        from = celltype_rename$old_name,
        to = celltype_rename$new_name
)
clin.assoc.pancurx$signature <- plyr::mapvalues(gsub(" Rep ", " ", clin.assoc.pancurx$signature),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
clin.assoc.pancurx <- clin.assoc.pancurx |>
        filter(!grepl("Ambient RNA|MALAT1|HSP|Proliferation", signature)) |>
        column_to_rownames("signature") |>
        mutate(log_HR_pancurx = `log(HR)`)

clin.assoc.pancurx.pval <- readRDS(snakemake@input[["clin_assoc_pancurx_pval"]]) |>
        mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
clin.assoc.pancurx.pval$signature <- plyr::mapvalues(gsub(" Rep ", " ", clin.assoc.pancurx.pval$signature),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
clin.assoc.pancurx.pval <- clin.assoc.pancurx.pval |>
        filter(!grepl("Ambient RNA|MALAT1|HSP|Proliferation", signature)) |>
        column_to_rownames("signature") |>
        mutate(log_HR_pancurx = `log(HR)`)


## niches
clin.assoc.tcga.niche <- readRDS(snakemake@input[["clin_assoc_paad_niche"]]) |>
        column_to_rownames("niche") |>
        mutate(celltype = NA) |>
        mutate(log_HR_tcga = `log(HR)`)
clin.assoc.pancurx.niche <- readRDS(snakemake@input[["clin_assoc_pancurx_niche"]]) |>
        column_to_rownames("niche") |>
        mutate(celltype = NA) |>
        mutate(log_HR_pancurx = `log(HR)`)

clin.assoc.tcga.pval.niche <- readRDS(snakemake@input[["clin_assoc_paad_pval_niche"]]) |>
        column_to_rownames("niche") |>
        mutate(log_HR_tcga = `log(HR)`)
clin.assoc.pancurx.pval.niche <- readRDS(snakemake@input[["clin_assoc_pancurx_pval_niche"]]) |>
        column_to_rownames("niche") |>
        mutate(log_HR_pancurx = `log(HR)`)

## construct matrices for plotting
mat1_sig <- cbind(clin.assoc.pancurx, clin.assoc.tcga) |>
        select(-celltype, -`log(HR)`) |>
        as.matrix()
mat1_niche <- cbind(clin.assoc.pancurx.niche, clin.assoc.tcga.niche) |>
        select(-celltype, -`log(HR)`) |>
        as.matrix()

mat2_sig <- cbind(clin.assoc.pancurx.pval, clin.assoc.tcga.pval) |>
        select(-`log(HR)`) |>
        as.matrix()
mat2_niche <- cbind(clin.assoc.pancurx.pval.niche, clin.assoc.tcga.pval.niche) |>
        select(-`log(HR)`) |>
        as.matrix()

column.split <- c(rep("Chan", length(clin.assoc.pancurx) - 3), rep("TCGA", length(clin.assoc.tcga) - 3))

mat1_combined <- rbind(
        mat1_sig[, !grepl("log_HR", colnames(mat1_sig))],
        mat1_niche[, !grepl("log_HR", colnames(mat1_niche))]
)
mat2_combined <- rbind(
        mat2_sig[, !grepl("log_HR", colnames(mat2_sig))],
        mat2_niche[, !grepl("log_HR", colnames(mat2_niche))]
)

mat1_combined_HR <- rbind(
        mat1_sig[, c("log_HR_pancurx", "log_HR_tcga")],
        mat1_niche[, c("log_HR_pancurx", "log_HR_tcga")]
)
mat2_combined_HR <- rbind(
        mat2_sig[, c("log_HR_pancurx", "log_HR_tcga")],
        mat2_niche[, c("log_HR_pancurx", "log_HR_tcga")]
)

row.split <- factor(c(clin.assoc.pancurx$celltype, rep("Ecosystem", nrow(mat1_niche)-1)), levels = c(unique(clin.assoc.pancurx$celltype), "Ecosystem"))

mat1_all <- cbind(mat1_combined, mat1_combined_HR)
mat2_all <- cbind(mat2_combined, mat2_combined_HR)

## plot the heatmap
colnames(mat1_all) <- plyr::mapvalues(colnames(mat1_all),
        from = c("log_HR_pancurx", "log_HR_tcga"),
        to = c("log(HR) - Chan", "log(HR) - TCGA")
)
colnames(mat1_all) <- gsub("treatment_outcome_first_course", "Treatment Outcome: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("AdjuvantBest_Therapy_Response", "Adjuvant Therapy Response: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("Complete Remission/Response", "CR", colnames(mat1_all))
colnames(mat1_all) <- gsub("Partial Remission/Response", "PR", colnames(mat1_all))
colnames(mat1_all) <- gsub("Stable Disease", "SD", colnames(mat1_all))
colnames(mat1_all) <- gsub("Progressive Disease", "PD", colnames(mat1_all))
colnames(mat1_all) <- gsub("Disease Progression", "PD", colnames(mat1_all))
colnames(mat1_all) <- gsub("new_tumor_event_type", "New Tumor Event: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("ajcc_pathologic_tumor_stage", "Pathological Stage: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("Clinical_Staging", "Clinical Stage: Stage ", colnames(mat1_all))
# colnames(mat1_all) <- gsub("PathologicalStaging", "Pathological Stage: Stage", colnames(mat1_all))
colnames(mat1_all) <- gsub("Histological_Grade", "Grade: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("histological_grade", "Grade: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("AdjuvantTreatmentType", "Adjuvant Therapy: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("Moffitt.mod.y", "Moffitt: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("basal-like", "Basal-like", colnames(mat1_all))
colnames(mat1_all) <- gsub("classic", "Classical", colnames(mat1_all))
colnames(mat1_all) <- gsub("PurIST", "PurIST: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("tumor_status", "Status: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("New Tumor Event", "New Event: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("new_tumor_event_site", "Metastasis: ", colnames(mat1_all))
colnames(mat1_all) <- gsub("ploidy", "Ploidy", colnames(mat1_all))
colnames(mat1_all) <- gsub("purity", "Purity", colnames(mat1_all))

rownames(mat1_all) <- plyr::mapvalues(rownames(mat1_all),
        from = c("Niche_1", "Niche_2", "Niche_3", "Niche_4"),
        to = c("Ecosystem0", "Ecosystem1 - Classical", "Ecosystem2 - Basal-like", "Ecosystem3 - Immune act.")
)

colnames(mat2_all) <- colnames(mat1_all)
rownames(mat2_all) <- rownames(mat1_all)

# remove "Ecosystem0" row
mat1_all <- mat1_all[!grepl("Ecosystem0", rownames(mat1_all)), ]
mat2_all <- mat2_all[!grepl("Ecosystem0", rownames(mat2_all)), ]

ht <- Heatmap(mat1_all,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        height = nrow(mat1_all) * unit(0.26, "in"),
        cell_fun = function(j, i, x, y, w, h, fill) {
                if (mat2_all[i, j] < 0.001) {
                        grid.text("***", x, y)
                } else if (mat2_all[i, j] < 0.01) {
                        grid.text("**", x, y)
                } else if (mat2_all[i, j] < 0.05) {
                        grid.text("*", x, y)
                }
        },
        #height = nrow(mat1_all) * unit(0.15, "in"),
        #width = ncol(mat1_all) * unit(0.21, "in"),
        name = "t-statistic/log(HR)",
        column_split = c(column.split, rep("", ncol(mat1_combined_HR))),
        row_split = row.split,
        cluster_row_slices = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 90,
        row_title = NULL,
        row_names_side = "left",
        row_dend_side = "right",
        left_annotation = rowAnnotation(
                Celltype = c(clin.assoc.pancurx$celltype, clin.assoc.pancurx.niche$celltype[1:3]), # tmp fix, removing Ecosystem0
                col = list(Celltype = celltype_pal_to_use),
                show_annotation_name = FALSE
        )
)

png(file = snakemake@output[["figure4_e_supp"]], width = 15, height = 20, units = "in", res = 600)
draw(ht)
dev.off()

margin_spacer <- function(x) {
  # where x is the column in your dataset
  left_length <- nchar(levels(factor(x)))[1]
  if (left_length > 8) {
    return((left_length - 8) * 4)
  }
  else
    return(0)
}

# Panel # - REMOVED (bar plots of signature/program HRs)
res.cox.to.plot <- readRDS(snakemake@input[["coxph_summary_paad"]]) |> 
        mutate(celltype = str_split(signature, "\\.", simplify = T)[, 1]) |>
        mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
res.cox.to.plot$celltype <- plyr::mapvalues(res.cox.to.plot$celltype,
        from = celltype_rename$old_name,
        to = celltype_rename$new_name
)
res.cox.to.plot$signature <- plyr::mapvalues(gsub(" Rep ", " ", res.cox.to.plot$signature),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
res.cox.to.plot <- res.cox.to.plot |>
        filter(!grepl("Ambient RNA|MALAT1|HSP|Proliferation", signature)) |>
        group_by(signature) |>
        mutate(unique_signature = paste0(signature, "_-_", row_number())) |>
        ungroup() |>
        arrange(desc(HR))
hr.paad <- ggplot(res.cox.to.plot, aes(x = reorder(unique_signature, -HR), y = HR)) +
        geom_bar(stat = "identity", aes(fill = celltype)) +
        scale_x_discrete(labels = res.cox.to.plot$signature) +
        geom_hline(yintercept = 1, linetype = "dashed") +
        scale_fill_manual(values = celltype_pal_to_use) +
        labs(x = NULL, y = "Hazard Ratio") +
        geom_text(aes(label = ifelse(p.adj < 0.05, "*", "")),
                position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt
        ) +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")
ggsave(snakemake@output[["figure4_bar_supp_tcga"]], width = 10, height = 8, units = "in", dpi = 600)
rm(res.cox.to.plot)

# Panel # - REMOVED (bar plots of signature/program HRs)
res.cox.to.plot <- readRDS(snakemake@input[["coxph_summary_pancurx"]]) |>
        mutate(celltype = str_split(signature, "\\.", simplify = T)[, 1]) |>
        mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
res.cox.to.plot$celltype <- plyr::mapvalues(res.cox.to.plot$celltype,
        from = celltype_rename$old_name,
        to = celltype_rename$new_name
)
res.cox.to.plot$signature <- plyr::mapvalues(gsub(" Rep ", " ", res.cox.to.plot$signature),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
res.cox.to.plot <- res.cox.to.plot |>
        filter(!grepl("Ambient RNA|MALAT1|HSP|Proliferation", signature)) |>
        group_by(signature) |>
        mutate(unique_signature = paste0(signature, "_-_", row_number())) |>
        ungroup() |>
        arrange(desc(HR))
hr.pancurx <- ggplot(res.cox.to.plot, aes(x = reorder(unique_signature, -HR), y = HR)) +
        geom_bar(stat = "identity", aes(fill = celltype)) +
        scale_x_discrete(labels = res.cox.to.plot$signature) +
        geom_hline(yintercept = 1, linetype = "dashed") +
        scale_fill_manual(values = celltype_pal_to_use) +
        labs(x = NULL, y = "Hazard Ratio") +
        geom_text(aes(label = ifelse(p.adj < 0.05, "*", "")),
                position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt
        ) +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")
ggsave(snakemake@output[["figure4_bar_supp_pancurx"]], width = 10, height = 8, units = "in", dpi = 600)
rm(res.cox.to.plot)

# Panel A - OLD, go to supp
print("Drawing Niche vs Epithelial subtype Program scatter plot")
paad.clin.for.plotting <- readRDS(snakemake@input[["clin_data_for_plotting_paad_niche"]])
pancurx.clin.for.plotting <- readRDS(snakemake@input[["clin_data_for_plotting_pancurx_niche"]])

paad.clin.for.plotting$Group <- "TCGA"
pancurx.clin.for.plotting$Group <- "Chan"

clin.for.plotting <- rbind(
        paad.clin.for.plotting |> 
          select(matches("^Group$") | contains("Niche_2_mean") | contains("Niche_3_mean") | contains("pancreatic epithelial cell.pancreatic epithelial cell Rep 4_mean") | contains("pancreatic epithelial cell.pancreatic epithelial cell Rep 14_mean")) #|>
          #pivot_longer(contains("Niche"), names_to = "Ecotype", values_to = "Ecotype score") |>
          #pivot_longer(contains("pancreatic epithelial cell"), names_to = "Epithelial Program", values_to = "Epithelial Program score")
          ,
        pancurx.clin.for.plotting |> 
          select(matches("^Group$") | contains("Niche_2_mean") | contains("Niche_3_mean") | contains("pancreatic epithelial cell.pancreatic epithelial cell Rep 4_mean") | contains("pancreatic epithelial cell.pancreatic epithelial cell Rep 14_mean")) #|>
          #pivot_longer(contains("Niche"), names_to = "Ecotype", values_to = "Ecotype score") |>
          #pivot_longer(contains("pancreatic epithelial cell"), names_to = "Epithelial Program", values_to = "Epithelial Program score")
        )


niche.2.vs.epi.14 <- ggplot(clin.for.plotting, aes(x = Niche_2_mean, y = `pancreatic epithelial cell.pancreatic epithelial cell Rep 14_mean`, colour = Group)) +
        geom_point() +
        scale_colour_jama() +
        lims(x = c(-30, NA)) +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor(show.legend = FALSE, method = "spearman", cor.coef.name = "rho") +
        labs(x = "Enrichment score - Ecosystem - Classical", y = "Enrichment score - Epithelial Classical-A", colour = "Cohort") +
        theme_pubr() +
        theme(legend.position = "top")

niche.3.vs.epi.4 <- ggplot(clin.for.plotting, aes(x = Niche_3_mean, y = `pancreatic epithelial cell.pancreatic epithelial cell Rep 4_mean`, colour = Group)) +
        geom_point() +
        scale_colour_jama() +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor(show.legend = FALSE, method = "spearman", cor.coef.name = "rho") +
        labs(x = "Enrichment score - Ecosystem - Basal-like", y = "Enrichment score - Epithelial Basal-like-A - EMT", colour = "Cohort") +
        theme_pubr() +
        theme(legend.position = "top")

niche.vs.epi <- ggarrange(niche.2.vs.epi.14, niche.3.vs.epi.4, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave(snakemake@output[["figure4_a_supp"]], plot = niche.vs.epi, width = snakemake@params[["figure4_a_width"]], height = snakemake@params[["figure4_a_height"]], units = "in", dpi = 600)
print("Niche vs Epithelial subtype Program scatter plot successfully created")

# Panel A
panel_a <- image_read_pdf(snakemake@input[['panel_a']]) |> image_ggplot()
ggsave(snakemake@output[["figure4_a"]], plot = panel_a, width = snakemake@params[["figure4_a_width"]], height = snakemake@params[["figure4_a_height"]], units = "in", dpi = 600)
print("KRAS variants vs. niche scores plot successfully created")

# Panel B - supp
panel_b <- image_read_pdf(snakemake@input[['panel_b']]) |> image_ggplot()
ggsave(snakemake@output[["figure4_b_supp"]], plot = panel_b, width = 7, height = 7, units = "in", dpi = 600)
print("Niche Loadings scatter plot successfully created")

# Panel B
print("Drawing Epithelial Program 10 survival curve")
p_list <- readRDS(snakemake@input[["cm_curve_plot_list_paad"]])
epi_10_cm_curve_tcga <- p_list[[19]]$plot + labs(title = "TCGA", color = "KRAS amp. program", fill = "KRAS amp. program")

p_list <- readRDS(snakemake@input[["cm_curve_plot_list_pancurx"]])
epi_10_cm_curve_pancurx <- p_list[[19]]$plot + labs(title = "Chan", color = "KRAS amp. program", fill = "KRAS amp. program")

epi_10_cm_curve <- ggarrange(epi_10_cm_curve_pancurx, epi_10_cm_curve_tcga, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")
ggsave(snakemake@output[["figure4_b"]], plot = ggarrange(epi_10_cm_curve_tcga, epi_10_cm_curve_pancurx, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right"), 
       width = snakemake@params[["figure4_b_width"]], height = snakemake@params[["figure4_b_height"]], units = "in", dpi = 600, bg = "white")
print("Epithelial Program 10 survival curve successfully created")

# Panel C
print("Drawing Chan Basal-like Niche vs Epithelial Program 10 scatter plot")
niche.3.vs.epi.10 <- ggplot(pancurx.clin.for.plotting, aes(x = Niche_3_mean, y = `pancreatic epithelial cell.pancreatic epithelial cell Rep 10_mean`)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor(method = "spearman", cor.coef.name = "rho") +
        labs(x = "Ecosystem - Basal-like", y = "Epithelial KRAS amp.", title = "Chan") +
        theme_pubr() +
        theme(legend.position = "top", plot.title = element_text(size = 18))
ggsave(snakemake@output[["figure4_c"]], plot = niche.3.vs.epi.10, width = snakemake@params[["figure4_c_width"]], height = snakemake@params[["figure4_c_height"]], units = "in", dpi = 600)
print("Chan Basal-like Niche vs Epithelial Program 10 scatter plot successfully created")

# Panel D
print("Drawing TCGA molecular subtype survival curve")
p_list <- readRDS(snakemake@input[["cm_curve_plot_list_paad_niche"]])
niche_3_cm_curve_tcga <- p_list[[6]]$plot

## new code for Panel D
niche_3_cm_curve_tcga <- survminer::ggsurvplot(
	fit = survival::survfit(survival::Surv(OS.time, OS) ~ PurIST, data = paad.clin.for.plotting), 
	data = paad.clin.for.plotting, 
	pval = TRUE, 
	conf.int = TRUE,
	legend.labs = c("Basal-like", "Classical"),  # Custom legend labels for TCGA
	palette = c("#f8766d", "#04bb36")
)
niche_3_cm_curve_tcga <- niche_3_cm_curve_tcga$plot + labs(title = "TCGA", color = "Subtype", fill = "Subtype") + theme(legend.position = "none")

#names(pancurx.clin.for.plotting)

niche_3_cm_curve_pancurx <- survminer::ggsurvplot(
	fit = survival::survfit(survival::Surv(OSFromSurgery, patient.vital_status) ~ `sscsubtype_coarse`, data = pancurx.clin.for.plotting), 
	data = pancurx.clin.for.plotting, 
	pval = TRUE, 
	conf.int = TRUE,
	legend.labs = c("Basal-like", "Classical", "Hybrid"),  # Custom legend labels for PanCurx
	palette = c("#f8766d", "#04bb36", "#619cff")
)
niche_3_cm_curve_pancurx <- niche_3_cm_curve_pancurx$plot + labs(title = "Chan", color = "Subtype", fill = "Subtype") + theme(legend.position = "none")

# p_list <- readRDS(snakemake@input[["cm_curve_plot_list_pancurx_niche"]])
# niche_3_cm_curve_pancurx <- p_list[[6]]$plot + 
# 	labs(title = "Chan\nExpression subtype") + 
# 	theme(title = element_text(size = 20)) +
# 	scale_color_discrete(labels = c("Basal-like", "Classical"))  # Custom legend labels for PanCurx

niche_3_cm_curve <- ggarrange(niche_3_cm_curve_pancurx, niche_3_cm_curve_tcga, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")

png(file = snakemake@output[["figure4_d"]], width = snakemake@params[["figure4_d_width"]], height = snakemake@params[["figure4_d_height"]], units = "in", res = 600)
niche_3_cm_curve
dev.off()
print("TCGA molecular subtype survival curve successfully created")

# Panel E - Chan
# Prepare data
col14_data <- data.frame(
  name = rownames(mat1_all),
  value = mat1_all[,14],
  pvalue = mat2_all[,14]
)

# Filter out Ecosystem0
col14_data <- col14_data[col14_data$name != "Ecosystem0", ]

# Add category and significance
col14_data$category <- ifelse(grepl("Ecosystem", col14_data$name), "Ecosystem", "Gene Program")
col14_data$stars <- ifelse(col14_data$pvalue < 0.001, "***",
                          ifelse(col14_data$pvalue < 0.01, "**",
                                 ifelse(col14_data$pvalue < 0.05, "*", "")))

# Order within categories
col14_data <- do.call(rbind, lapply(split(col14_data, col14_data$category), function(x) {
  x$name <- factor(x$name, levels = x$name[order(x$value, decreasing = TRUE)])
  return(x)
}))

write_tsv(col14_data, snakemake@output[["figure4_e_chan_data"]])

# Create plot
panel_e_chan <- ggplot(col14_data, aes(x = name, y = -log10(pvalue))) +
  geom_bar(stat = "identity", 
           aes(fill = interaction(category, name == "Ecosystem2 - Basal-like"))) +
  scale_fill_manual(values = c("grey60", "red", "lightblue", "lightblue")) +
  #geom_text(aes(label = stars), vjust = -0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_grid(~category, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "-log10(p-value)",
    title = "Chan cohort - Adjuvant Therapy - PD"
  )
ggsave(snakemake@output[["figure4_e_chan"]], plot = panel_e_chan, width = snakemake@params[["figure4_e_width"]], height = snakemake@params[["figure4_e_height"]], units = "in", dpi = 600, bg = "white")

panel_e_chan <- ggplot(col14_data, aes(x = category, y = value)) +
  geom_hline(yintercept=0, linetype=3, colour='grey60') +
  geom_violindot(data = filter(col14_data, category == "Gene Program"), fill='#c3d9db') +
  geom_point(data = filter(col14_data, category == "Ecosystem"), color='#8f7d7c', size=5) +
  theme_pubr() +
  labs(y = expression("Association with therapy response (worse outcome " %->% ")"),
       subtitle = "Chan (Tumour LCM enriched)") +
  theme(axis.title.x = element_blank()) +
  geom_text(data = filter(col14_data, category == "Ecosystem"), aes(label = name), hjust = 0.4, vjust = 2) #+
#   coord_flip() +
#   theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
#         axis.title.y = element_blank()) +
#   labs(x = "Association with therapy response (worse outcome →)")
  

# Panel E - TCGA
# Prepare data
col34_data <- data.frame(
  name = rownames(mat1_all),
  value = mat1_all[,34],
  pvalue = mat2_all[,34]
)

# Filter out Ecosystem0
col34_data <- col34_data[col34_data$name != "Ecosystem0", ]

# Add category and significance
col34_data$category <- ifelse(grepl("Ecosystem", col34_data$name), "Ecosystem", "Gene Program")
col34_data$stars <- ifelse(col34_data$pvalue < 0.001, "***",
			  ifelse(col34_data$pvalue < 0.01, "**",
				 ifelse(col34_data$pvalue < 0.05, "*", "")))

# Order within categories
col34_data <- do.call(rbind, lapply(split(col34_data, col34_data$category), function(x) {
  x$name <- factor(x$name, levels = x$name[order(x$value, decreasing = TRUE)])
  return(x)
}))

write_tsv(col34_data, snakemake@output[["figure4_e_tcga_data"]])

# Create plot
panel_e_tcga <- ggplot(col34_data, aes(x = name, y = -log10(pvalue))) +
  geom_bar(stat = "identity", 
	   aes(fill = interaction(category, name == "Ecosystem2 - Basal-like"))) +
  scale_fill_manual(values = c("grey60", "red", "lightblue", "lightblue")) +
  #geom_text(aes(label = stars), vjust = -0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_grid(~category, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "-log10(p-value)",
    title = "TCGA cohort - Therapy outcome - PD"
  )
ggsave(snakemake@output[["figure4_e_tcga"]], plot = panel_e_tcga, width = snakemake@params[["figure4_e_width"]], height = snakemake@params[["figure4_e_height"]], units = "in", dpi = 600, bg = "white")

panel_e_tcga <- ggplot(col34_data, aes(x = category, y = value)) +
  geom_hline(yintercept=0, linetype=3, colour='grey60') +
  geom_violindot(data = filter(col34_data, category == "Gene Program"), fill='#c3d9db') +
  geom_point(data = filter(col34_data, category == "Ecosystem"), color='#8f7d7c', size=5) +
  theme_pubr() +
  labs(y = expression("Association with therapy response (worse outcome " %->% ")"),
       subtitle = "TCGA (non tumour enriched)") +
  theme(axis.title.x = element_blank()) +
  geom_text(data = filter(col34_data, category == "Ecosystem"), aes(label = name), , hjust = 0.4, vjust = -1.5) #+
#   coord_flip() +
#   theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
#   	axis.title.y = element_blank()) +
#   labs(x = "Association with therapy response (better outcome →)")

# Panel F
# Prepare data for treatment outcomes
treatment_data <- data.frame(
  outcome = c("CR", "PR", "SD", "PD"),
  coefficient = mat1_all["Ecosystem3 - Immune act.", 
                        c("Treatment Outcome: CR", 
                          "Treatment Outcome: PR",
                          "Treatment Outcome: SD", 
                          "Treatment Outcome: PD")],
  pvalue = mat2_all["Ecosystem3 - Immune act.", 
                    c("Treatment Outcome: CR", 
                      "Treatment Outcome: PR",
                      "Treatment Outcome: SD", 
                      "Treatment Outcome: PD")]
)
treatment_data$outcome <- factor(treatment_data$outcome, 
                               levels = c("CR", "PR", "SD", "PD"))

# Add significance stars
treatment_data$stars <- ifelse(treatment_data$pvalue < 0.001, "***",
                              ifelse(treatment_data$pvalue < 0.01, "**",
                                    ifelse(treatment_data$pvalue < 0.05, "*", "")))

print(treatment_data)

# Create plot
panel_f <- ggplot(treatment_data, aes(x = outcome, y = coefficient)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = stars), hjust = -0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Treatment Outcome",
    y = "t-statistic",
    title = "TCGA\nEcosystem3 - Immune act."
  ) + 
  coord_flip()

# Panel G
# Prepare data for treatment outcomes
treatment_data <- data.frame(
  outcome = c("NED", "PD"),
  coefficient = mat1_all["Ecosystem3 - Immune act.", 
			c("Adjuvant Therapy Response: NED", 
			  "Adjuvant Therapy Response: PD")],
  pvalue = mat2_all["Ecosystem3 - Immune act.", 
		    c("Adjuvant Therapy Response: NED", 
		      "Adjuvant Therapy Response: PD")]
)
treatment_data$outcome <- factor(treatment_data$outcome, 
			       levels = c("NED", "PD"))

# Add significance stars
treatment_data$stars <- ifelse(treatment_data$pvalue < 0.001, "***",
			      ifelse(treatment_data$pvalue < 0.01, "**",
				    ifelse(treatment_data$pvalue < 0.05, "*", "")))
print(treatment_data)

# Create plot
panel_g <- ggplot(treatment_data, aes(x = outcome, y = coefficient)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = stars), hjust = -0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Adj. Therapy Response",
    y = "t-statistic",
    title = "Chan\nEcosystem3 - Immune act."
  ) + 
  coord_flip()

# Panel H
# Prepare data for treatment outcomes
treatment_data <- data.frame(
  outcome = c("CR", "PR", "SD", "PD"),
  coefficient = mat1_all["Ecosystem2 - Basal-like", 
			c("Treatment Outcome: CR", 
			  "Treatment Outcome: PR",
			  "Treatment Outcome: SD", 
			  "Treatment Outcome: PD")],
  pvalue = mat2_all["Ecosystem2 - Basal-like", 
		    c("Treatment Outcome: CR", 
		      "Treatment Outcome: PR",
		      "Treatment Outcome: SD", 
		      "Treatment Outcome: PD")]
)
treatment_data$outcome <- factor(treatment_data$outcome, 
			       levels = c("CR", "PR", "SD", "PD"))

# Add significance stars
treatment_data$stars <- ifelse(treatment_data$pvalue < 0.001, "***",
			      ifelse(treatment_data$pvalue < 0.01, "**",
				    ifelse(treatment_data$pvalue < 0.05, "*", "")))

print(treatment_data)

# Create plot
panel_h <- ggplot(treatment_data, aes(x = outcome, y = coefficient)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = stars), hjust = -0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Treatment Outcome",
    y = "t-statistic",
    title = "TCGA\nEcosystem2 - Basal-like"
  ) + 
  coord_flip()

# Panel I
# Prepare data for treatment outcomes
treatment_data <- data.frame(
  outcome = c("NED", "PD"),
  coefficient = mat1_all["Ecosystem2 - Basal-like", 
			c("Adjuvant Therapy Response: NED", 
			  "Adjuvant Therapy Response: PD")],
  pvalue = mat2_all["Ecosystem2 - Basal-like", 
		    c("Adjuvant Therapy Response: NED", 
		      "Adjuvant Therapy Response: PD")]
)
treatment_data$outcome <- factor(treatment_data$outcome, 
			       levels = c("NED", "PD"))

# Add significance stars
treatment_data$stars <- ifelse(treatment_data$pvalue < 0.001, "***",
			      ifelse(treatment_data$pvalue < 0.01, "**",
				    ifelse(treatment_data$pvalue < 0.05, "*", "")))

print(treatment_data)

# Create plot
panel_i <- ggplot(treatment_data, aes(x = outcome, y = coefficient)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = stars), hjust = -0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Adj. Therapy Response",
    y = "t-statistic",
    title = "Chan\nEcosystem2 - Basal-like"
  ) + 
  coord_flip()

# plot Figure 4

ht.grob <- grid.grabExpr(draw(ht))

design <- "
AAAABBBB
AAAABBBB
CCDDEEEE
CCDDEEEE
CCDDEEEE
CCDDEEEE
"

pdf(file = snakemake@output[["figure4_pdf"]], width = 13, height = 9)
epi_10_cm_curve + niche_3_cm_curve + panel_e_chan + panel_e_tcga + panel_a + 
        plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()

print("Figure 4 NEW PDF successfully created")

design <- "
AAAABB
AAAABB
AAAABB
AAAABB
AAAABB
CCCCDD
CCCCDD
CCCCDD
CCCCDD
CCCCDD
EEIIJJ
EEIIJJ
FFIIJJ
GGIIJJ
GGIIJJ
HHIIJJ
"

png(file = snakemake@output[["figure4_png"]], width = snakemake@params[["figure4_width"]], height = snakemake@params[["figure4_height"]], units = "in", res = 600)
epi_10_cm_curve + niche.3.vs.epi.10 + niche_3_cm_curve + panel_a + panel_h + panel_i + panel_f + panel_g + panel_e_chan + panel_e_tcga + 
        plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()

print("Figure 4 NEW PNG successfully created")