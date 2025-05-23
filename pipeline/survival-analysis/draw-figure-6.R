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

# Panel E
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
# names(clin.assoc.pancurx) <- gsub("^grade$", "Histological_Grade", names(clin.assoc.pancurx))
# names(clin.assoc.tcga) <- gsub("^grade$", "histological_grade", names(clin.assoc.tcga))
# names(clin.assoc.pancurx.niche) <- gsub("^grade$", "Histological_Grade", names(clin.assoc.pancurx.niche))
# names(clin.assoc.tcga.niche) <- gsub("^grade$", "histological_grade", names(clin.assoc.tcga.niche))

# names(clin.assoc.pancurx.pval) <- gsub("^grade$", "Histological_Grade", names(clin.assoc.pancurx.pval))
# names(clin.assoc.tcga.pval) <- gsub("^grade$", "histological_grade", names(clin.assoc.tcga.pval))
# names(clin.assoc.pancurx.pval.niche) <- gsub("^grade$", "Histological_Grade", names(clin.assoc.pancurx.pval.niche))
# names(clin.assoc.tcga.pval.niche) <- gsub("^grade$", "histological_grade", names(clin.assoc.tcga.pval.niche))

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
        scale(mat1_sig[, !grepl("log_HR", colnames(mat1_sig))]),
        scale(mat1_niche[, !grepl("log_HR", colnames(mat1_niche))])
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
        name = "log(HR)",
        cluster_columns = F,
        column_names_rot = 45,
        row_names_side = "left",
        row_dend_side = "right"
)

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
        name = "lm coefficient/log(HR)",
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

png(file = snakemake@output[["figure6_e"]], width = snakemake@params[["figure6_e_width"]], height = snakemake@params[["figure6_e_height"]], units = "in", res = 600)
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
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
ggsave("output/v3/figures/survival-analysis/figure-6/test-tcga-coxph-HR.png", width = 15, height = 13, units = "in", dpi = 600)
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
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
ggsave("output/v3/figures/survival-analysis/figure-6/test-pancurx-coxph-HR.png", width = 15, height = 13, units = "in", dpi = 600)
rm(res.cox.to.plot)

# Panel B
print("Drawing Epithelial Program 10 survival curve")
p_list <- readRDS(snakemake@input[["cm_curve_plot_list_paad"]])
epi_10_cm_curve_tcga <- p_list[[19]]$plot + labs(title = "TCGA - Epithelial Program\nZFAS1/P4HA1/EIF4A2")

p_list <- readRDS(snakemake@input[["cm_curve_plot_list_pancurx"]])
epi_10_cm_curve_pancurx <- p_list[[19]]$plot + labs(title = "Chan - Epithelial Program\nZFAS1/P4HA1/EIF4A2")

epi_10_cm_curve <- ggarrange(epi_10_cm_curve_pancurx, epi_10_cm_curve_tcga, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave(snakemake@output[["figure6_b"]], plot = ggarrange(epi_10_cm_curve_tcga, epi_10_cm_curve_pancurx, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right"), 
       width = snakemake@params[["figure6_b_width"]], height = snakemake@params[["figure6_b_height"]], units = "in", dpi = 600, bg = "white")
print("Epithelial Program 10 survival curve successfully created")

# Panel A
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
ggsave(snakemake@output[["figure6_a"]], plot = niche.vs.epi, width = snakemake@params[["figure6_a_width"]], height = snakemake@params[["figure6_a_height"]], units = "in", dpi = 600)
print("Niche vs Epithelial subtype Program scatter plot successfully created")

# Panel C
print("Drawing Chan Basal-like Niche vs Epithelial Program 10 scatter plot")
niche.3.vs.epi.10 <- ggplot(pancurx.clin.for.plotting, aes(x = Niche_3_mean, y = `pancreatic epithelial cell.pancreatic epithelial cell Rep 10_mean`)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor(method = "spearman", cor.coef.name = "rho") +
        labs(x = "Ecosystem - Basal-like", y = "Epithelial Program ZFAS1/P4HA1/EIF4A2", title = "Chan\nEnrichment score") +
        theme_pubr() +
        theme(legend.position = "top", plot.title = element_text(size = 18))
ggsave(snakemake@output[["figure6_c"]], plot = niche.3.vs.epi.10, width = snakemake@params[["figure6_c_width"]], height = snakemake@params[["figure6_c_height"]], units = "in", dpi = 600)
print("Chan Basal-like Niche vs Epithelial Program 10 scatter plot successfully created")

# Panel D
print("Drawing TCGA molecular subtype survival curve")
p_list <- readRDS(snakemake@input[["cm_curve_plot_list_paad_niche"]])
niche_3_cm_curve_tcga <- p_list[[6]]$plot

## new code for Panel D
niche_3_cm_curve_tcga <- survminer::ggsurvplot(fit = survival::survfit(survival::Surv(OS.time, OS) ~ PurIST, data = paad.clin.for.plotting), data = paad.clin.for.plotting, pval = TRUE, conf.int = TRUE, pval.coord = c(1800, 0.95))
niche_3_cm_curve_tcga <- niche_3_cm_curve_tcga$plot + labs(title = "TCGA\nExpression subtype")

p_list <- readRDS(snakemake@input[["cm_curve_plot_list_pancurx_niche"]])
niche_3_cm_curve_pancurx <- p_list[[6]]$plot + labs(title = "Chan\nNiche - Basal-like") + theme(title = element_text(size = 20))

png(file = snakemake@output[["figure6_d"]], width = snakemake@params[["figure6_d_width"]], height = snakemake@params[["figure6_d_height"]], units = "in", res = 600)
ggarrange(niche_3_cm_curve_tcga, niche_3_cm_curve_pancurx, ncol = 2, nrow = 1, common.legend = TRUE)
dev.off()
print("TCGA molecular subtype survival curve successfully created")

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