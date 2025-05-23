# fix some weird issue with the R environment
#BiocManager::install(version = "3.18")
.libPaths()
# Add custom library path if needed
#.libPaths(c(.libPaths(), "/ddn_exa/campbell/cyu/cyu/R/x86_64-pc-linux-gnu-library/4.2"))
#library(matrixStats, lib.loc = "/ddn_exa/campbell/cyu/cyu/R/x86_64-pc-linux-gnu-library/4.2")
sessionInfo()

#install.packages("survminer", repos = "https://cloud.r-project.org")
#remotes::install_version("matrixStats", version = "1.1.0", repos = "https://cloud.r-project.org", lib = "/home/campbell/cyu/R/x86_64-pc-linux-gnu-library/4.2")
#install.packages("Rcpp", repos = "https://cloud.r-project.org", INSTALL_opts = '--no-lock')
#library(Rcpp)

# Load libraries
suppressPackageStartupMessages({
        library(tidyverse)
        library(magrittr)
        library(sjstats)
        library(reticulate)
        library(fastcluster)
        library(BiocParallel)
        library(GSVA)
        library(survival)
        library(survminer)
        library(readxl)
})

score_for_clin_assoc <- snakemake@params[["score_used_for_clinical_association"]]
score_for_clin_assoc <- paste0("_", score_for_clin_assoc)

# load signature interpt and cell type rename
sig.interpt <- read_xlsx(snakemake@input[["sig_interpretation"]])
celltype_rename <- read_csv(snakemake@input[["cell_type_rename"]])

# load PDAC signatures
genesets.sig <- lapply(snakemake@input[["celltype_gene_loading_mtxs"]], read_tsv)

names(genesets.sig) <- snakemake@params[["celltypes"]]

genesets.sig <- lapply(genesets.sig, function(df) {
        w <- df %>%
                distinct(gene, .keep_all = T) %>%
                column_to_rownames("gene")

        l <- lapply(w, function(sig) {
                ll <- rownames(w)[sort(sig, index.return = TRUE, decreasing = TRUE)$ix[1:50]]
                str_split(ll, "_ENSG", simplify = T)[, 1]
        })
        names(l) <- names(w)
        l
})
genesets.sig <- unlist(genesets.sig, recursive = FALSE, use.names = TRUE)


# TCGA survival analysis
## load TCGA expression and clinical data
PAAD.mRNA.vst <- readRDS(snakemake@input[["tidy_paad_mrna_data"]])
paad.clin <- read_tsv(snakemake@input[["tidy_paad_clin_data"]])

print("Data loaded")

## GSVA for the signatures in TCGA
#unloadNamespace("matrixStats")
#remotes::install_version("matrixStats", version = "1.1.0")
#library(matrixStats)

expr_matrix <- PAAD.mRNA.vst |> column_to_rownames("bcr_patient_barcode") |> as.matrix() |> t() # t() is needed to make it genes x samples
# Create GSVA params object
gsva_params <- gsvaParam(
    exprData = expr_matrix,
    geneSets = genesets.sig,
    kcdf = "Poisson"
)
# Run GSVA
gsva.es <- gsva(
        param = gsva_params,
        verbose = TRUE,
        BPPARAM = MulticoreParam(16)
) |> t() |> as.data.frame() |> rownames_to_column("bcr_patient_barcode") # t() is needed to make it samples x genes
# gsva.es <- gsva(
#         expr = PAAD.mRNA.vst %>% column_to_rownames("bcr_patient_barcode") %>% as.matrix() %>% t(),
#         gset.idx.list = genesets.sig,
#         method = "gsva",
#         kcdf = "Poisson",
#         verbose = T,
#         BPPARAM = MulticoreParam(16)
# ) %>%
#         t() %>%
#         as.data.frame() %>%
#         rownames_to_column("bcr_patient_barcode")
print("GSVA done")

## survival analysis
paad.clin.for.plotting <- paad.clin
p_list <- list()
clin.assoc.list <- list()
cox.list <- list()
cox.gsva.list <- list()
res.cut.list <- list()
res.cut.gsva.list <- list()

for (gs in names(genesets.sig)) {
        # Take the mRNA data
        PAAD.mRNA.vst.geneset <- PAAD.mRNA.vst %>%
                # then make it a tibble (nice printing while debugging)
                as_tibble() %>%
                # then get just a few genes
                select(bcr_patient_barcode, all_of(intersect(names(.), genesets.sig[[gs]]))) %>%
                # then trim the barcode (see head(clin), and ?substr)
                mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
                # then z-scale the expressions
                mutate(across(where(is.numeric), scale)) %>%
                # then sum/mean over genes
                mutate(
                        sum = rowSums(across(where(is.numeric))),
                        mean = rowMeans(across(where(is.numeric)))
                ) %>%
                # then get GSVA/ssGSEA enrichment scores
                mutate(gsva = gsva.es[[gs]]) %>%
                # get high/low exprs
                mutate(
                        geneset_status = ifelse(mean > median(mean), "Top 50%", "Bottom 50%"),
                        geneset_status_gsva = ifelse(gsva > median(gsva), "high", "low")
                ) %>%
                # subset to wanted columns
                select(bcr_patient_barcode, sum, mean, gsva, geneset_status, geneset_status_gsva)

        # then join back to clinical data
        paad.clin.geneset <- inner_join(PAAD.mRNA.vst.geneset, paad.clin, by = "bcr_patient_barcode")

        names(PAAD.mRNA.vst.geneset) <- c("bcr_patient_barcode", paste0(gs, "_sum"), paste0(gs, "_mean"), paste0(gs, "_gsva"), paste0(gs, "_status"), paste0(gs, "_status_gsva"))
        paad.clin.for.plotting <- inner_join(PAAD.mRNA.vst.geneset, paad.clin.for.plotting, by = "bcr_patient_barcode")

        # further remove samples with missing/not useful information
        paad.clin.geneset <- paad.clin.geneset[!is.na(paad.clin.geneset$PFI.time), ]
        #print("#############Survival analysis for geneset#############")
        #print(gs)
        #print(head(paad.clin.geneset |> select(bcr_patient_barcode, PFI.time, PFI, mean, gsva, geneset_status, geneset_status_gsva)))

        # model association with clinical variables
        # clin.assoc.list[["gender"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + gender")), data = paad.clin.for.plotting))$coefficients
        clin.assoc.list[["grade"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + histological_grade")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("G4|GX", histological_grade)) #|>
                        # mutate(histological_grade = plyr::mapvalues(histological_grade, 
                        #                                             from = c("G1", "G2", "G3"), 
                        #                                             to = c("1", "2", "3")) |> as.numeric())
        ))$coefficients
        clin.assoc.list[["stage"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + ajcc_pathologic_tumor_stage")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("Discrepancy|^Stage I$", ajcc_pathologic_tumor_stage))
        ))$coefficients
        clin.assoc.list[["status"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + tumor_status")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("Discrepancy", tumor_status))
        ))$coefficients
        clin.assoc.list[["newtumor"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + new_tumor_event_type")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("Primary|\\|", new_tumor_event_type))
        ))$coefficients
        clin.assoc.list[["newtumorsite"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + new_tumor_event_site")),
                data = paad.clin.for.plotting |>
                        filter(grepl("Liver|Lung", new_tumor_event_site))
        ))$coefficients
        clin.assoc.list[["outcome"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + treatment_outcome_first_course")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("Discrepancy", treatment_outcome_first_course))
        ))$coefficients
        clin.assoc.list[["subtype"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + PurIST")), data = paad.clin.for.plotting))$coefficients
        # clin.assoc.list[["subtype2"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + Moffitt")), data = paad.clin.for.plotting))$coefficients

        # fit CoxPH model
        cox1 <- coxph(Surv(PFI.time, PFI) ~ mean, data = paad.clin.geneset)
        cox2 <- coxph(Surv(PFI.time, PFI) ~ gsva, data = paad.clin.geneset)
        cox.list[[gs]] <- cox1
        cox.gsva.list[[gs]] <- cox2

        # find optimal cutpoints to categorize geneset loading
        res.cut <- surv_cutpoint(paad.clin.geneset,
                time = "PFI.time", event = "PFI",
                variables = c("mean"), minprop = 0.4
        )
        res.cut.gsva <- surv_cutpoint(paad.clin.geneset,
                time = "PFI.time", event = "PFI",
                variables = c("gsva"), minprop = 0.4
        )
        #print("#############Cutpoint#############")
        #summary(res.cut) |> print()
        #print("#############Cutpoint GSVA#############")
        #summary(res.cut.gsva) |> print()
        res.cut.list[[gs]] <- res.cut
        res.cut.gsva.list[[gs]] <- res.cut.gsva

        # plot(res.cut, "mean", palette = "npg")
        # plot(res.cut, "gsva", palette = "npg")

        res.cat <- surv_categorize(res.cut)
        res.cat.gsva <- surv_categorize(res.cut.gsva)
        #print("#############Survival categorize#############")
        #head(res.cat) |> print()
        #print("#############Survival categorize GSVA#############")
        #head(res.cat.gsva) |> print()

        # create survival curves
        sfit1 <- survfit(Surv(PFI.time, PFI) ~ mean, data = res.cat)
        sfit2 <- survfit(Surv(PFI.time, PFI) ~ gsva, data = res.cat.gsva)
        #sfit1 <- survfit(Surv(PFI.time, PFI) ~ geneset_status, data = paad.clin.geneset)
        #sfit2 <- survfit(Surv(PFI.time, PFI) ~ geneset_status_gsva, data = paad.clin.geneset)
        #print("#############Survival fit#############")
        #summary(sfit1, times=seq(0,365*5,365)) |> print()
        #print("#############Survival fit GSVA#############")
        #summary(sfit2, times=seq(0,365*5,365)) |> print()

        celltype <- plyr::mapvalues(str_split(gs, "\\.", simplify = TRUE)[, 1],
                from = celltype_rename$old_name,
                to = celltype_rename$new_name,
                warn_missing = FALSE
        )
        sig <- plyr::mapvalues(gsub(" Rep ", " ", str_split(gs, "\\.", simplify = TRUE)[, 2]),
                from = sig.interpt$signature,
                to = sig.interpt$`short interpretation`,
                warn_missing = FALSE
        )

        p_list <- c(
                p_list,
                ggsurvplot_list(
                        fit = list(sfit1, sfit2),
                        #data = paad.clin.geneset,
                        data = list(res.cat, res.cat.gsva),
                        legend.title = list("Expression quantile", "GSVA"),
                        conf.int = TRUE, pval = TRUE, risk.table = TRUE, title = paste0(celltype, " - ", sig),
                        #pval.coord = c(1800, .95), # p-value location
                        palette = c("#E7B800", "#2E9FDF")
                )
        )
}
rm(gs, celltype, sig)

# save survival curves
png("output/v3/figures/survival-analysis/test.png", width = 10, height = 10, units = "in", res = 600)
p_list[[19]]$plot + labs(x = "Time (days)", y = "Survival probability") + theme(legend.position = "none")
dev.off()

print("Survival analysis done")

## Extract data from cox.list
univ_results <- lapply(cox.list, function(x) {
        x <- summary(x)
        p.value <- x$waldtest["pvalue"]
        wald.test <- x$waldtest["test"]
        beta <- x$coefficients[1] # coefficient beta
        HR <- x$coefficients[2] # exp(beta)
        HR.confint.lower <- x$conf.int[, "lower .95"]
        HR.confint.upper <- x$conf.int[, "upper .95"]
        HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
        res <- c(beta, HR, wald.test, p.value)
        names(res) <- c(
                "beta", "HR (95% CI for HR)", "wald.test",
                "p.value"
        )
        return(res)
})
res.cox <- t(as.data.frame(univ_results, check.names = FALSE))
res.cox <- as.data.frame(res.cox)

res.cox <- res.cox |>
        mutate(p.value = as.numeric(p.value)) |>
        mutate(p.adj = p.adjust(p.value, method = "fdr")) |>
        mutate(HR = str_split(res.cox$`HR (95% CI for HR)`, " ", simplify = T)[, 1] |> as.numeric()) |>
        rownames_to_column("signature") #|>
        #mutate(celltype = str_split(signature, "\\.", simplify = T)[, 1]) |>
        #mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])

# res.cox$celltype <- plyr::mapvalues(res.cox$celltype,
#         from = celltype_rename$old_name,
#         to = celltype_rename$new_name
# )
# res.cox$signature <- plyr::mapvalues(gsub(" Rep ", " ", res.cox$signature),
#         from = sig.interpt$signature,
#         to = sig.interpt$`short interpretation`
# )

res.cox.tcga <- res.cox

## Extract data from clin.assoc.list
holder.list <- lapply(clin.assoc.list, function(clin.assoc) {
        clin.assoc <- lapply(clin.assoc, function(x) x[, 3]) # 1 for estimated coefficient, 3 for t-statistic
        as.data.frame(do.call(rbind, clin.assoc)) |> rownames_to_column("signature")
})
holder.list <- lapply(names(holder.list), function(clin.var) {
        df <- holder.list[[clin.var]]
 	#df <- df |> mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
 	# df$signature <- plyr::mapvalues(gsub(" Rep ", " ", df$signature),
        #  	from = sig.interpt$signature,
        #  	to = sig.interpt$`short interpretation`
 	# )
        names(df) <- plyr::mapvalues(names(df), from = c("V1"), to = c(clin.var))
        df
})
names(holder.list) <- names(clin.assoc.list)

clin.assoc.tcga <- Reduce(function(x, y) {
        left_join(x, y, by = "signature")
}, holder.list)
holder <- res.cox |>
        dplyr::select(signature, HR) |>
        mutate(HR = log(HR))
names(holder) <- c("signature", "log(HR)")
clin.assoc.tcga <- left_join(holder, clin.assoc.tcga, by = "signature")

rm(holder.list, holder)

holder.list <- lapply(clin.assoc.list, function(clin.assoc) {
        clin.assoc <- lapply(clin.assoc, function(x) x[, 4])
        as.data.frame(do.call(rbind, clin.assoc)) %>% rownames_to_column("signature")
})
holder.list <- lapply(names(holder.list), function(clin.var) {
        df <- holder.list[[clin.var]]
 	#df <- df |> mutate(signature = str_split(signature, "\\.", simplify = T)[, 2])
 	# df$signature <- plyr::mapvalues(gsub(" Rep ", " ", df$signature),
        #  	from = sig.interpt$signature,
        #  	to = sig.interpt$`short interpretation`
 	# )
        names(df) <- plyr::mapvalues(names(df), from = c("V1"), to = c(clin.var))
        df
})
names(holder.list) <- names(clin.assoc.list)

clin.assoc.tcga.pval <- Reduce(function(x, y) {
        left_join(x, y, by = "signature")
}, holder.list)
holder <- res.cox |> select(signature, p.adj)
names(holder) <- c("signature", "log(HR)")
clin.assoc.tcga.pval <- left_join(holder, clin.assoc.tcga.pval, by = "signature")

rm(holder.list, holder)


# TO SAVE
saveRDS(gsva.es, snakemake@output[["gsva_es_paad"]])
saveRDS(paad.clin.for.plotting, snakemake@output[["clin_data_for_plotting_paad"]])
saveRDS(p_list, snakemake@output[["cm_curve_plot_list_paad"]])
saveRDS(res.cox.tcga, snakemake@output[["coxph_summary_paad"]])
saveRDS(clin.assoc.tcga, snakemake@output[["clin_assoc_paad"]])
saveRDS(clin.assoc.tcga.pval, snakemake@output[["clin_assoc_paad_pval"]])


# PanCuRx survival analysis
## load PanCuRx expression and clinical data
pancurx.mRNA.logTPM <- readRDS(snakemake@input[["tidy_pancurx_mrna_data"]])
pancurx.clin.new <- read_tsv(snakemake@input[["tidy_pancurx_clin_data"]])

print("Data loaded")

## GSVA for the signatures in PanCuRx
# unloadNamespace("matrixStats")
# remotes::install_version("matrixStats", version = "1.1.0")
# library(matrixStats)

expr_matrix <- pancurx.mRNA.logTPM |> column_to_rownames("sample_id") |> as.matrix() |> t() # t() is needed to make it genes x samples
# Create GSVA params object
gsva_params <- gsvaParam(
    exprData = expr_matrix,
    geneSets = genesets.sig,
    kcdf = "Poisson"
)
# Run GSVA
gsva.es <- gsva(
        param = gsva_params,
        verbose = TRUE,
        BPPARAM = MulticoreParam(16)
) |> t() |> as.data.frame() |> rownames_to_column("sample_id") # t() is needed to make it samples x genes

# gsva.es <- gsva(
#         expr = pancurx.mRNA.logTPM %>% column_to_rownames("sample_id") %>% as.matrix() %>% t(),
#         gset.idx.list = genesets.sig,
#         method = "gsva",
#         kcdf = "Poisson",
#         verbose = T,
#         BPPARAM = MulticoreParam(16)
# ) %>%
#         t() %>%
#         as.data.frame() %>%
#         rownames_to_column("sample_id")
print("GSVA done")

## survival analysis
pancurx.clin.for.plotting <- pancurx.clin.new
p_list <- list()
clin.assoc.list <- list()
cox.list <- list()
cox.gsva.list <- list()
res.cut.list <- list()
res.cut.gsva.list <- list()

for (gs in names(genesets.sig)) {
        # Take the mRNA data
        pancurx.mRNA.logTPM.geneset <- pancurx.mRNA.logTPM %>%
                # then make it a tibble (nice printing while debugging)
                as_tibble() %>%
                # then get just a few genes
                select(sample_id, all_of(intersect(names(.), genesets.sig[[gs]]))) %>%
                # then trim the barcode (see head(clin), and ?substr)
                # mutate(sample_id = substr(sample_id, 1, 12)) %>%
                # then z-scale the expressions
                mutate(across(where(is.numeric), scale)) %>%
                # then sum/mean over genes
                mutate(
                        sum = rowSums(across(where(is.numeric))),
                        mean = rowMeans(across(where(is.numeric)))
                ) %>%
                # then get GSVA/ssGSEA enrichment scores
                mutate(gsva = gsva.es[[gs]]) %>%
                # get high/low exprs
                mutate(
                        geneset_status = ifelse(mean > median(mean), "Top 50%", "Bottom 50%"),
                        geneset_status_gsva = ifelse(gsva > median(gsva), "high", "low")
                ) %>%
                # subset to wanted columns
                select(sample_id, sum, mean, gsva, geneset_status, geneset_status_gsva)

        # then join back to clinical data
        pancurx.clin.geneset <- inner_join(pancurx.mRNA.logTPM.geneset, pancurx.clin.new, by = "sample_id")

        names(pancurx.mRNA.logTPM.geneset) <- c("sample_id", paste0(gs, "_sum"), paste0(gs, "_mean"), paste0(gs, "_gsva"), paste0(gs, "_status"), paste0(gs, "_status_gsva"))
        pancurx.clin.for.plotting <- inner_join(pancurx.mRNA.logTPM.geneset, pancurx.clin.for.plotting, by = "sample_id")

        # further remove samples with missing/not useful information
        pancurx.clin.geneset <- pancurx.clin.geneset[!is.na(pancurx.clin.geneset$OSFromSurgery), ]
        # pancurx.clin.geneset <- pancurx.clin.geneset[!grepl("Unknown", pancurx.clin.geneset$AdjuvantTreatmentType),]

        # model association with clinical variables
        clin.assoc.list[["ploidy"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + ploidy")), data = pancurx.clin.for.plotting))$coefficients
        clin.assoc.list[["purity"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + purity")), data = pancurx.clin.for.plotting))$coefficients
        # clin.assoc.list[["sex"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + Sex")), data = pancurx.clin.for.plotting))$coefficients
        clin.assoc.list[["stage"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + Clinical_Staging")),
                data = pancurx.clin.for.plotting |>
                        #filter(!grepl("^I$|III", Clinical_Staging))
                        #mutate(Clinical_Staging = gsub("A$|B$", "", Clinical_Staging))
                        filter(!grepl("^I$|^II$", Clinical_Staging))
        ))$coefficients
        # clin.assoc.list[["stage2"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + PathologicalStaging")),
        #         data = pancurx.clin.for.plotting |>
        #                 filter(!grepl("I - III", PathologicalStaging)) |>
        #                 mutate(PathologicalStaging = gsub("^IA$|^IB$", "I", PathologicalStaging)) |>
        #                 mutate(PathologicalStaging = gsub("^III$|^IV$", "III - IV", PathologicalStaging))
        # ))$coefficients
        clin.assoc.list[["grade"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + Histological_Grade")),
                data = pancurx.clin.for.plotting |>
                        filter(!grepl("N/A|Other|Unknown", Histological_Grade)) |>
                        mutate(Histological_Grade = plyr::mapvalues(Histological_Grade, 
                                                                    from = c("Well differentiated", "Moderately differentiated", "Poorly differentiated", "Undifferentiated"), 
                                                                    to = c("G1", "G2", "G3", "G4")))
        ))$coefficients
        clin.assoc.list[["adjuvant"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + AdjuvantTreatmentType")),
                data = pancurx.clin.for.plotting |>
                        filter(!grepl("Unknown", AdjuvantTreatmentType))
        ))$coefficients
        clin.assoc.list[["adjuvant_outcome"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + AdjuvantBest_Therapy_Response")),
                data = pancurx.clin.for.plotting |>
                        filter(!grepl("N/A|Unknown", AdjuvantBest_Therapy_Response))
        ))$coefficients
        # clin.assoc.list[["macrophage"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + CD68.PPC.Pos.Pix.Perc.Stroma")),
        #         data = pancurx.clin.for.plotting |>
        #                 filter(!grepl("Discrepancy", CD68.PPC.Pos.Pix.Perc.Stroma))
        # ))$coefficients
        clin.assoc.list[["subtype"]][[gs]] <- summary(lm(as.formula(paste0("`", gs, score_for_clin_assoc, "` ~ 0 + Moffitt.mod.y")),
                data = pancurx.clin.for.plotting
        ))$coefficients

        # fit CoxPH model
        cox1 <- coxph(Surv(OSFromSurgery, patient.vital_status) ~ mean, data = pancurx.clin.geneset)
        cox2 <- coxph(Surv(OSFromSurgery, patient.vital_status) ~ gsva, data = pancurx.clin.geneset)
        # cox1 <- coxph(Surv(OSFromSurgery, patient.vital_status) ~ mean + purity + Sex + AdjuvantTreatmentType, data = pancurx.clin.geneset)
        # cox2 <- coxph(Surv(OSFromSurgery, patient.vital_status) ~ gsva + purity + Sex + AdjuvantTreatmentType, data = pancurx.clin.geneset)
        cox.list[[gs]] <- cox1
        cox.gsva.list[[gs]] <- cox2

        # find optimal cutpoints to categorize geneset loading
        res.cut <- surv_cutpoint(pancurx.clin.geneset,
                time = "OSFromSurgery", event = "patient.vital_status",
                variables = c("mean"), minprop = 0.4
        )
        res.cut.gsva <- surv_cutpoint(pancurx.clin.geneset,
                time = "OSFromSurgery", event = "patient.vital_status",
                variables = c("gsva"), minprop = 0.4
        )
        # summary(res.cut)
        res.cut.list[[gs]] <- res.cut
        res.cut.gsva.list[[gs]] <- res.cut.gsva

        # plot(res.cut, "mean", palette = "npg")
        # plot(res.cut, "gsva", palette = "npg")

        res.cat <- surv_categorize(res.cut)
        res.cat.gsva <- surv_categorize(res.cut.gsva)
        # head(res.cat)

        # create survival curves
        # sfit1 <- survfit(Surv(OSFromSurgery, patient.vital_status) ~ geneset_status, data = pancurx.clin.geneset)
        # sfit2 <- survfit(Surv(OSFromSurgery, patient.vital_status) ~ geneset_status_gsva, data = pancurx.clin.geneset)
        sfit1 <- survfit(Surv(OSFromSurgery, patient.vital_status) ~ mean, data = res.cat)
        sfit2 <- survfit(Surv(OSFromSurgery, patient.vital_status) ~ gsva, data = res.cat.gsva)
        # summary(sfit, times=seq(0,365*5,365))

        celltype <- plyr::mapvalues(str_split(gs, "\\.", simplify = TRUE)[, 1],
                from = celltype_rename$old_name,
                to = celltype_rename$new_name,
                warn_missing = FALSE
        )
        sig <- plyr::mapvalues(gsub(" Rep ", " ", str_split(gs, "\\.", simplify = TRUE)[, 2]),
                from = sig.interpt$signature,
                to = sig.interpt$`short interpretation`,
                warn_missing = FALSE
        )

        p_list <- c(
                p_list,
                ggsurvplot_list(list(sfit1, sfit2),
                        # data = pancurx.clin.geneset,
                        data = list(res.cat, res.cat.gsva),
                        legend.title = list("Expression quantile", "gsva"),
                        conf.int = TRUE, pval = TRUE, risk.table = TRUE, title = paste0(celltype, " - ", sig),
                        #pval.coord = c(2800, .95), # p-value location
                        palette = c("#E7B800", "#2E9FDF")
                )
        )
}
rm(gs, celltype, sig)

print("Survival analysis done")

## Extract data from cox.list
univ_results <- lapply(cox.list, function(x) {
        x <- summary(x)
        p.value <- x$waldtest["pvalue"]
        wald.test <- x$waldtest["test"]
        beta <- x$coefficients[1] # coefficient beta
        HR <- x$coefficients[2] # exp(beta)
        HR.confint.lower <- x$conf.int[, "lower .95"]
        HR.confint.upper <- x$conf.int[, "upper .95"]
        HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
        res <- c(beta, HR, wald.test, p.value)
        names(res) <- c(
                "beta", "HR (95% CI for HR)", "wald.test",
                "p.value"
        )
        return(res)
})
res.cox <- t(as.data.frame(univ_results, check.names = FALSE))
res.cox <- as.data.frame(res.cox)
rm(univ_results)

res.cox <- res.cox |>
        mutate(p.value = as.numeric(p.value)) |>
        mutate(p.adj = p.adjust(p.value, method = "fdr")) |>
        mutate(HR = str_split(res.cox$`HR (95% CI for HR)`, " ", simplify = T)[, 1] |> as.numeric()) |>
        rownames_to_column("signature") #|>
	# mutate(celltype = str_split(signature, "\\.", simplify = T)[,1]) |>
	# mutate(signature = str_split(signature, "\\.", simplify = T)[,2])

# res.cox$celltype <- plyr::mapvalues(res.cox$celltype,
#         from = celltype_rename$old_name,
#         to = celltype_rename$new_name
# )
# res.cox$signature <- plyr::mapvalues(gsub(" Rep ", " ", res.cox$signature),
#         from = sig.interpt$signature,
#         to = sig.interpt$`short interpretation`
# )

res.cox.pancurx <- res.cox
print("#############CoxPH summary#############")
print(res.cox.pancurx)

## Extract data from clin.assoc.list
holder.list <- lapply(clin.assoc.list, function(clin.assoc) {
        clin.assoc <- lapply(clin.assoc, function(x) x[, 3]) # 1 for estimated coefficient, 3 for t-statistic
        as.data.frame(do.call(rbind, clin.assoc)) %>% rownames_to_column("signature")
})
# holder.list <- lapply(names(holder.list), function(clin.var) {
#         df <- holder.list[[clin.var]]
#         names(df) <- plyr::mapvalues(names(df), from = c("V1"), to = c(clin.var))
#         df
# })
names(holder.list) <- names(clin.assoc.list)

clin.assoc.pancurx <- Reduce(function(x, y) {
        left_join(x, y, by = "signature")
}, holder.list)
holder <- res.cox |>
        select(signature, HR) |>
        mutate(HR = log(HR))
names(holder) <- c("signature", "log(HR)")
clin.assoc.pancurx <- left_join(holder, clin.assoc.pancurx, by = "signature")
rm(holder.list, holder)

print("#############Clinical association#############")
print(clin.assoc.pancurx)

holder.list <- lapply(clin.assoc.list, function(clin.assoc) {
        clin.assoc <- lapply(clin.assoc, function(x) x[, 4])
        as.data.frame(do.call(rbind, clin.assoc)) %>% rownames_to_column("signature")
})
holder.list <- lapply(names(holder.list), function(clin.var) {
        df <- holder.list[[clin.var]]
        names(df) <- plyr::mapvalues(names(df), from = c("V1"), to = c(clin.var))
        df
})
names(holder.list) <- names(clin.assoc.list)

clin.assoc.pancurx.pval <- Reduce(function(x, y) {
        left_join(x, y, by = "signature")
}, holder.list)
holder <- res.cox |> select(signature, p.adj)
names(holder) <- c("signature", "log(HR)")
clin.assoc.pancurx.pval <- left_join(holder, clin.assoc.pancurx.pval, by = "signature")
rm(holder.list, holder)

print("#############Clinical association p-value#############")
print(clin.assoc.pancurx.pval)

# TO SAVE
saveRDS(gsva.es, snakemake@output[["gsva_es_pancurx"]])
saveRDS(pancurx.clin.for.plotting, snakemake@output[["clin_data_for_plotting_pancurx"]])
saveRDS(p_list, snakemake@output[["cm_curve_plot_list_pancurx"]])
saveRDS(res.cox.pancurx, snakemake@output[["coxph_summary_pancurx"]])
saveRDS(clin.assoc.pancurx, snakemake@output[["clin_assoc_pancurx"]])
saveRDS(clin.assoc.pancurx.pval, snakemake@output[["clin_assoc_pancurx_pval"]])

