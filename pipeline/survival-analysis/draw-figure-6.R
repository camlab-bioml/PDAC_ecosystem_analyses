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

# Panel A
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
        mutate(HR_minus_1_tcga = `HR-1`)

clin.assoc.tcga.pval <- readRDS(snakemake@input[["clin_assoc_paad_pval"]]) |> 
        mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
clin.assoc.tcga.pval$signature <- plyr::mapvalues(gsub(" Rep ", " ", clin.assoc.tcga.pval$signature),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
clin.assoc.tcga.pval <- clin.assoc.tcga.pval |>
        filter(!grepl("Ambient RNA|MALAT1|HSP|Proliferation", signature)) |>
        column_to_rownames("signature") |>
        mutate(HR_minus_1_tcga = `HR-1`)

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
        mutate(HR_minus_1_pancurx = `HR-1`)

clin.assoc.pancurx.pval <- readRDS(snakemake@input[["clin_assoc_pancurx_pval"]]) |>
        mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
clin.assoc.pancurx.pval$signature <- plyr::mapvalues(gsub(" Rep ", " ", clin.assoc.pancurx.pval$signature),
        from = sig.interpt$signature,
        to = sig.interpt$`short interpretation`
)
clin.assoc.pancurx.pval <- clin.assoc.pancurx.pval |>
        filter(!grepl("Ambient RNA|MALAT1|HSP|Proliferation", signature)) |>
        column_to_rownames("signature") |>
        mutate(HR_minus_1_pancurx = `HR-1`)


## niches
clin.assoc.tcga.niche <- readRDS(snakemake@input[["clin_assoc_paad_niche"]]) |>
        column_to_rownames("niche") |>
        mutate(celltype = NA) |>
        mutate(HR_minus_1_tcga = `HR-1`)
clin.assoc.pancurx.niche <- readRDS(snakemake@input[["clin_assoc_pancurx_niche"]]) |>
        column_to_rownames("niche") |>
        mutate(celltype = NA) |>
        mutate(HR_minus_1_pancurx = `HR-1`)

clin.assoc.tcga.pval.niche <- readRDS(snakemake@input[["clin_assoc_paad_pval_niche"]]) |>
        column_to_rownames("niche") |>
        mutate(HR_minus_1_tcga = `HR-1`)
clin.assoc.pancurx.pval.niche <- readRDS(snakemake@input[["clin_assoc_pancurx_pval_niche"]]) |>
        column_to_rownames("niche") |>
        mutate(HR_minus_1_pancurx = `HR-1`)

## construct matrices for plotting
# names(clin.assoc.pancurx) <- gsub("^grade$", "Histological_Grade", names(clin.assoc.pancurx))
# names(clin.assoc.tcga) <- gsub("^grade$", "histological_grade", names(clin.assoc.tcga))
# names(clin.assoc.pancurx.niche) <- gsub("^grade$", "Histological_Grade", names(clin.assoc.pancurx.niche))
# names(clin.assoc.tcga.niche) <- gsub("^grade$", "histological_grade", names(clin.assoc.tcga.niche))

# names(clin.assoc.pancurx.pval) <- gsub("^grade$", "Histological_Grade", names(clin.assoc.pancurx.pval))
# names(clin.assoc.tcga.pval) <- gsub("^grade$", "histological_grade", names(clin.assoc.tcga.pval))
# names(clin.assoc.pancurx.pval.niche) <- gsub("^grade$", "Histological_Grade", names(clin.assoc.pancurx.pval.niche))
# names(clin.assoc.tcga.pval.niche) <- gsub("^grade$", "histological_grade", names(clin.assoc.tcga.pval.niche))

mat1_sig <- cbind(clin.assoc.pancurx, clin.assoc.tcga) |>
        select(-celltype, -`HR-1`) |>
        as.matrix()
mat1_niche <- cbind(clin.assoc.pancurx.niche, clin.assoc.tcga.niche) |>
        select(-celltype, -`HR-1`) |>
        as.matrix()

mat2_sig <- cbind(clin.assoc.pancurx.pval, clin.assoc.tcga.pval) |>
        select(-`HR-1`) |>
        as.matrix()
mat2_niche <- cbind(clin.assoc.pancurx.pval.niche, clin.assoc.tcga.pval.niche) |>
        select(-`HR-1`) |>
        as.matrix()

column.split <- c(rep("Chan", length(clin.assoc.pancurx) - 3), rep("TCGA", length(clin.assoc.tcga) - 3))

mat1_combined <- rbind(
        scale(mat1_sig[, !grepl("HR_minus_1", colnames(mat1_sig))]),
        scale(mat1_niche[, !grepl("HR_minus_1", colnames(mat1_niche))])
)

mat2_combined <- rbind(
        mat2_sig[, !grepl("HR_minus_1", colnames(mat2_sig))],
        mat2_niche[, !grepl("HR_minus_1", colnames(mat2_niche))]
)

mat1_combined_HR <- rbind(
        mat1_sig[, c("HR_minus_1_pancurx", "HR_minus_1_tcga")],
        mat1_niche[, c("HR_minus_1_pancurx", "HR_minus_1_tcga")]
)

mat2_combined_HR <- rbind(
        mat2_sig[, c("HR_minus_1_pancurx", "HR_minus_1_tcga")],
        mat2_niche[, c("HR_minus_1_pancurx", "HR_minus_1_tcga")]
)

#row.split <- c(rep("Signature", nrow(mat1_sig)), rep("Niche", nrow(mat1_niche)))
row.split <- factor(c(clin.assoc.pancurx$celltype, rep("Ecotype", nrow(mat1_niche))), levels = c(unique(clin.assoc.pancurx$celltype), "Ecotype"))

mat1_all <- cbind(mat1_combined, mat1_combined_HR)
mat2_all <- cbind(mat2_combined, mat2_combined_HR)

## plot the heatmap
h1 <- Heatmap(mat1_combined,
        cell_fun = function(j, i, x, y, w, h, fill) {
                if (mat2_combined[i, j] < 0.001) {
                        grid.text("***", x, y)
                } else if (mat2_combined[i, j] < 0.01) {
                        grid.text("**", x, y)
                } else if (mat2_combined[i, j] < 0.05) {
                        grid.text("*", x, y)
                }
        },
        height = nrow(mat1_combined) * unit(0.15, "in"),
        width = ncol(mat1_combined) * unit(0.21, "in"),
        name = "lm coefficient",
        column_split = column.split,
        row_split = row.split,
        cluster_columns = F,
        column_names_rot = 45,
        row_title = NULL,
        row_names_side = "left",
        row_dend_side = "right",
        left_annotation = rowAnnotation(
                Celltype = c(clin.assoc.pancurx$celltype, clin.assoc.pancurx.niche$celltype),
                col = list(Celltype = celltype_pal_to_use),
                show_annotation_name = FALSE
        )
)
h2 <- Heatmap(mat1_combined_HR,
        cell_fun = function(j, i, x, y, w, h, fill) {
                if (mat2_combined_HR[i, j] < 0.001) {
                        grid.text("***", x, y)
                } else if (mat2_combined_HR[i, j] < 0.01) {
                        grid.text("**", x, y)
                } else if (mat2_combined_HR[i, j] < 0.05) {
                        grid.text("*", x, y)
                }
        },
        height = nrow(mat1_combined_HR) * unit(0.15, "in"),
        width = nrow(mat1_combined_HR) * unit(0.21, "in"),
        name = "HR - 1",
        cluster_columns = F,
        column_names_rot = 45,
        row_names_side = "left",
        row_dend_side = "right"
)

colnames(mat1_all) <- plyr::mapvalues(colnames(mat1_all),
        from = c("HR_minus_1_pancurx", "HR_minus_1_tcga"),
        to = c("HR - 1 (Chan)", "HR - 1 (TCGA)")
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
        to = c("Ecotype - Unknown", "Ecotype 1 - Classical", "Ecotype 2 - Basal-like", "Ecotype 3 - Immune act.")
)

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
        name = "lm coefficient/HR - 1",
        column_split = c(column.split, rep("", ncol(mat1_combined_HR))),
        row_split = row.split,
        cluster_row_slices = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 90,
        row_title = NULL,
        row_names_side = "left",
        row_dend_side = "right",
        left_annotation = rowAnnotation(
                Celltype = c(clin.assoc.pancurx$celltype, clin.assoc.pancurx.niche$celltype),
                col = list(Celltype = celltype_pal_to_use),
                show_annotation_name = FALSE
        )
)

margin_spacer <- function(x) {
  # where x is the column in your dataset
  left_length <- nchar(levels(factor(x)))[1]
  if (left_length > 8) {
    return((left_length - 8) * 4)
  }
  else
    return(0)
}

# Panel B
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
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
#ggsave("test-tcga-coxph-HR.png", width = 15, height = 13, units = "in", dpi = 321)
rm(res.cox.to.plot)

# Panel C
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
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
# ggsave("test-pancurx-coxph-HR.png", width = 15, height = 13, units = "in", dpi = 321)
rm(res.cox.to.plot)

# Panel D
p_list <- readRDS(snakemake@input[["cm_curve_plot_list_paad"]])

epi_10_cm_curve_tcga <- p_list[[19]]$plot + labs(title = "TCGA - Epithelial Program\nZFAS1/P4HA1/EIF4A2")

# Panel E
p_list <- readRDS(snakemake@input[["cm_curve_plot_list_pancurx"]])

epi_10_cm_curve_pancurx <- p_list[[19]]$plot + labs(title = "Chan - Epithelial Program\nZFAS1/P4HA1/EIF4A2")

# Panel F
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
        lims(x = c(-30, NA)) +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor(show.legend = FALSE) +
        labs(x = "Ecotype - Classical", y = "Epithelial Classical-A", title = "Enrichment score") +
        theme_pubr() +
        theme(legend.position = "top")

niche.3.vs.epi.4 <- ggplot(clin.for.plotting, aes(x = Niche_3_mean, y = `pancreatic epithelial cell.pancreatic epithelial cell Rep 4_mean`, colour = Group)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor() +
        labs(x = "Ecotype - Basal-like", y = "Epithelial Basal-like-A - EMT", title = "Enrichment score") +
        theme_pubr() +
        theme(legend.position = "top")

niche.vs.epi <- ggarrange(niche.2.vs.epi.14, niche.3.vs.epi.4, ncol = 2, nrow = 1, common.legend = TRUE)

niche.3.vs.epi.10 <- ggplot(pancurx.clin.for.plotting, aes(x = Niche_3_mean, y = `pancreatic epithelial cell.pancreatic epithelial cell Rep 10_mean`)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor() +
        labs(x = "Ecotype - Basal-like", y = "Epithelial Program ZFAS1/P4HA1/EIF4A2", title = "Chan\nEnrichment score") +
        theme_pubr() +
        theme(legend.position = "top")

# Panel G
p_list <- readRDS(snakemake@input[["cm_curve_plot_list_paad_niche"]])
niche_3_cm_curve_tcga <- p_list[[6]]$plot

## new code for panel D
niche_3_cm_curve_tcga <- survminer::ggsurvplot(fit = survival::survfit(survival::Surv(OS.time, OS) ~ PurIST, data = paad.clin.for.plotting), data = paad.clin.for.plotting, pval = TRUE) 
niche_3_cm_curve_tcga <- niche_3_cm_curve_tcga$plot + labs(title = "TCGA\nExpression subtype")

p_list <- readRDS(snakemake@input[["cm_curve_plot_list_pancurx_niche"]])
niche_3_cm_curve_pancurx <- p_list[[6]]$plot + labs(title = "Chan\nEcotype - Basal-like")

niche_3_cm_curve <- ggarrange(niche_3_cm_curve_tcga, niche_3_cm_curve_pancurx, ncol = 2, nrow = 1, common.legend = TRUE)

# plot Figure 2
design <- "
AAAAEEEEEE
AAAAEEEEEE
BBBBEEEEEE
BBBBEEEEEE
CCDDEEEEEE
CCDDEEEEEE
"

ht.grob <- grid.grabExpr(draw(ht))
epi_10_cm_curve <- ggarrange(epi_10_cm_curve_tcga, epi_10_cm_curve_pancurx, ncol = 2, nrow = 1, common.legend = TRUE)

pdf(file = snakemake@output[["figure6_pdf"]], width = snakemake@params[["figure6_width"]], height = snakemake@params[["figure6_height"]])
niche.vs.epi + epi_10_cm_curve + niche.3.vs.epi.10 + niche_3_cm_curve_tcga + ht.grob + #hr.paad + hr.pancurx + 
        plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()

print("Figure 6 PDF successfully created")

png(file = snakemake@output[["figure6_png"]], width = snakemake@params[["figure6_width"]], height = snakemake@params[["figure6_height"]], units = "in", res = 360)
niche.vs.epi + epi_10_cm_curve + niche.3.vs.epi.10 + niche_3_cm_curve_tcga + ht.grob + #hr.paad + hr.pancurx + 
        plot_layout(design = design) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()

print("Figure 6 PNG successfully created")