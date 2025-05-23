# Load libraries
suppressPackageStartupMessages({
        library(magrittr)
        library(tidyverse)
        library(sjstats)
        library(tidyr)
        library(stringr)
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

# load nich loadings in patients
niche.factors <- readRDS(snakemake@input[["niche_factor_loadings"]])

niche.factors <- t(scale(t(niche.factors)))
colnames(niche.factors) <- gsub(" ", "_", colnames(niche.factors))
rownames(niche.factors) <- gsub("  ", " ", rownames(niche.factors))

# TCGA survival analysis
## load TCGA clinical data for plotting
paad.clin.for.plotting <- readRDS(snakemake@input[["clin_data_for_plotting_paad"]])

## GSVA for the signatures in TCGA
gsva.es <- readRDS(snakemake@input[["gsva_es_paad"]])

## niche enrichment scores - GSVA
niche.gsva.es <- gsva.es |>
        column_to_rownames("bcr_patient_barcode")

names(niche.gsva.es) <- str_split(names(niche.gsva.es), "\\.", simplify = T)[, 2]
names(niche.gsva.es) <- gsub(" Rep", "", names(niche.gsva.es))
#names(niche.gsva.es) <- gsub("-|, ", "_", names(niche.gsva.es))
names(niche.gsva.es) <- gsub(" ", "_", names(niche.gsva.es))

print("Checking if the columns in niche.factors are present in niche.gsva.es")
table(colnames(niche.factors) %in% names(niche.gsva.es))
print("Checking if the columns in niche.gsva.es are present in niche.factors")
table(names(niche.gsva.es) %in% colnames(niche.factors))

niche.gsva.es <- niche.gsva.es[, colnames(niche.factors)]
niche.gsva.es <- as.matrix(niche.gsva.es)

niche.gsva.es <- niche.gsva.es %*% t(niche.factors)
niche.gsva.es <- as.data.frame(niche.gsva.es) |>
        rownames_to_column("bcr_patient_barcode")
names(niche.gsva.es) <- gsub(" ", "_", names(niche.gsva.es))

## niche enrichment scores - marker expression mean values
head(paad.clin.for.plotting)
niche.sig.exprs.mean <- paad.clin.for.plotting |>
        select(1, ends_with("_mean")) |>
        column_to_rownames("bcr_patient_barcode")

names(niche.sig.exprs.mean) <- str_split(names(niche.sig.exprs.mean), "\\.", simplify = T)[, 2]
names(niche.sig.exprs.mean) <- gsub(" Rep", "", names(niche.sig.exprs.mean))
#names(niche.sig.exprs.mean) <- gsub("-|, ", "_", names(niche.sig.exprs.mean))
names(niche.sig.exprs.mean) <- gsub(" ", "_", names(niche.sig.exprs.mean))
names(niche.sig.exprs.mean) <- gsub("_mean", "", names(niche.sig.exprs.mean))

print("Checking if the columns in niche.factors are present in niche.sig.exprs.mean")
table(colnames(niche.factors) %in% names(niche.sig.exprs.mean))
print("Checking if the columns in niche.sig.exprs.mean are present in niche.factors")
table(names(niche.sig.exprs.mean) %in% colnames(niche.factors))

niche.sig.exprs.mean <- niche.sig.exprs.mean[, colnames(niche.factors)]
niche.sig.exprs.mean <- as.matrix(niche.sig.exprs.mean)

niche.sig.exprs.mean <- niche.sig.exprs.mean %*% t(niche.factors)
niche.sig.exprs.mean <- as.data.frame(niche.sig.exprs.mean) |>
        rownames_to_column("bcr_patient_barcode")
names(niche.sig.exprs.mean) <- gsub(" ", "_", names(niche.sig.exprs.mean))


## survival analysis
p_list <- list()
clin.assoc.list <- list()
cox.list <- list()
cox.gsva.list <- list()
res.cut.list <- list()
res.cut.gsva.list <- list()

for (niche in grep("Niche", colnames(niche.gsva.es), value = T)) {
        niche.gsva.es.subset <- niche.gsva.es |>
                mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) |>
                select(c("bcr_patient_barcode", niche))
        names(niche.gsva.es.subset) <- c("bcr_patient_barcode", "gsva")

	niche.sig.exprs.mean.subset <- niche.sig.exprs.mean |>
		mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) |>
		select(c("bcr_patient_barcode", niche))
	names(niche.sig.exprs.mean.subset) <- c("bcr_patient_barcode", "mean")

        # join niche gsva estimate and expression score to clinical data
        paad.clin.geneset <- inner_join(niche.gsva.es.subset, paad.clin.for.plotting, by = "bcr_patient_barcode")
	paad.clin.geneset <- inner_join(niche.sig.exprs.mean.subset, paad.clin.geneset, by = "bcr_patient_barcode")

        names(niche.gsva.es.subset) <- c("bcr_patient_barcode", paste0(niche, "_gsva"))
	names(niche.sig.exprs.mean.subset) <- c("bcr_patient_barcode", paste0(niche, "_mean"))
        paad.clin.for.plotting <- inner_join(niche.gsva.es.subset, paad.clin.for.plotting, by = "bcr_patient_barcode")
	paad.clin.for.plotting <- inner_join(niche.sig.exprs.mean.subset, paad.clin.for.plotting, by = "bcr_patient_barcode")

        # further remove samples with missing/not useful information
        paad.clin.geneset <- paad.clin.geneset[!is.na(paad.clin.geneset$PFI.time), ]

        # model association with clinical variables
        # clin.assoc.list[["gender"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + gender")), data = paad.clin.for.plotting))$coefficients
        clin.assoc.list[["grade"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + histological_grade")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("G4|GX", histological_grade)) #|>
                        # mutate(histological_grade = plyr::mapvalues(histological_grade, 
                        #                                             from = c("G1", "G2", "G3"), 
                        #                                             to = c("1", "2", "3")) |> as.numeric())
        ))$coefficients
        clin.assoc.list[["stage"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + ajcc_pathologic_tumor_stage")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("Discrepancy|^Stage I$", ajcc_pathologic_tumor_stage))
        ))$coefficients
        clin.assoc.list[["status"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + tumor_status")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("Discrepancy", tumor_status))
        ))$coefficients
        clin.assoc.list[["newtumor"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + new_tumor_event_type")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("Primary|\\|", new_tumor_event_type))
        ))$coefficients
        clin.assoc.list[["newtumorsite"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + new_tumor_event_site")),
                data = paad.clin.for.plotting |>
                        filter(grepl("Liver|Lung", new_tumor_event_site))
        ))$coefficients
        clin.assoc.list[["outcome"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + treatment_outcome_first_course")),
                data = paad.clin.for.plotting |>
                        filter(!grepl("Discrepancy", treatment_outcome_first_course))
        ))$coefficients
        clin.assoc.list[["subtype"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + PurIST")), data = paad.clin.for.plotting))$coefficients
        # clin.assoc.list[["subtype2"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + Moffitt")), data = paad.clin.for.plotting))$coefficients

        # fit CoxPH model
  	cox1 <- coxph(Surv(PFI.time, PFI) ~ mean, data = paad.clin.geneset)
  	cox2 <- coxph(Surv(PFI.time, PFI) ~ gsva, data = paad.clin.geneset)
 	cox.list[[niche]] <- cox1
 	cox.gsva.list[[niche]] <- cox2

        # find optimal cutpoints to categorize geneset loading
 	res.cut <- surv_cutpoint(paad.clin.geneset,
         	time = "PFI.time", event = "PFI",
         	variables = c("mean"), minprop = 0.5
 	)
        res.cut.gsva <- surv_cutpoint(paad.clin.geneset,
                time = "PFI.time", event = "PFI",
                variables = c("gsva"), minprop = 0.5
        )
        # summary(res.cut)
 	res.cut.list[[niche]] <- res.cut
        res.cut.gsva.list[[niche]] <- res.cut.gsva

        # plot(res.cut.gsva, "niche_gsva", palette = "npg")
 	
	res.cat <- surv_categorize(res.cut)
        res.cat.gsva <- surv_categorize(res.cut.gsva)
        # head(res.cat.gsva)

        # create survival curves
 	sfit1 <- survfit(Surv(PFI.time, PFI) ~ mean, data = res.cat)
 	sfit2 <- survfit(Surv(PFI.time, PFI) ~ gsva, data = res.cat.gsva)
        # summary(sfit, times=seq(0,365*5,365))

        p_list <- c(
                p_list,
                ggsurvplot_list(list(sfit1, sfit2),
                        data = list(res.cat, res.cat.gsva),
                        legend.title = list("exprs", "gsva"),
                        conf.int = TRUE, pval = TRUE, risk.table = TRUE, title = niche
                )
        )
}
rm(niche)

## Extract data from cox.list
univ_results <- lapply(cox.list, function(x) {
        x <- summary(x)
        p.value <- x$waldtest["pvalue"]
        wald.test <- x$waldtest["test"]
        beta <- x$coefficients[1, 1] # coefficient beta
        HR <- x$coefficients[1, 2] # exp(beta)
        HR.confint.lower <- x$conf.int[1, "lower .95"]
        HR.confint.upper <- x$conf.int[1, "upper .95"]
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
        rownames_to_column("niche")

res.cox.tcga.niche <- res.cox

## Extract data from clin.assoc.list
holder.list <- lapply(clin.assoc.list, function(clin.assoc) {
        clin.assoc <- lapply(clin.assoc, function(x) x[, 3]) # 1 for estimate, 2 for std.error, 3 for t value, 4 for p value
        as.data.frame(do.call(rbind, clin.assoc)) %>% rownames_to_column("niche")
})
holder.list <- lapply(names(holder.list), function(clin.var) {
        df <- holder.list[[clin.var]]
        names(df) <- plyr::mapvalues(names(df), from = c("V1"), to = c(clin.var))
        df
})
names(holder.list) <- names(clin.assoc.list)

clin.assoc.tcga.niche <- Reduce(function(x, y) {
        left_join(x, y, by = "niche")
}, holder.list)
holder <- res.cox |>
        select(niche, HR) |>
        mutate(HR = log(HR))
names(holder) <- c("niche", "log(HR)")
clin.assoc.tcga.niche <- left_join(holder, clin.assoc.tcga.niche, by = "niche")

rm(holder.list, holder)

holder.list <- lapply(clin.assoc.list, function(clin.assoc) {
        clin.assoc <- lapply(clin.assoc, function(x) x[, 4])
        as.data.frame(do.call(rbind, clin.assoc)) %>% rownames_to_column("niche")
})
holder.list <- lapply(names(holder.list), function(clin.var) {
        df <- holder.list[[clin.var]]
        names(df) <- plyr::mapvalues(names(df), from = c("V1"), to = c(clin.var))
        df
})
names(holder.list) <- names(clin.assoc.list)

clin.assoc.tcga.pval.niche <- Reduce(function(x, y) {
        left_join(x, y, by = "niche")
}, holder.list)
holder <- res.cox |> select(niche, p.adj)
names(holder) <- c("niche", "log(HR)")
clin.assoc.tcga.pval.niche <- left_join(holder, clin.assoc.tcga.pval.niche, by = "niche")

rm(holder.list, holder)


# TO SAVE
saveRDS(niche.gsva.es, snakemake@output[["gsva_es_paad_niche"]])
saveRDS(niche.sig.exprs.mean, snakemake@output[["sig_exprs_mean_paad_niche"]])
saveRDS(paad.clin.for.plotting, snakemake@output[["clin_data_for_plotting_paad_niche"]])
saveRDS(p_list, snakemake@output[["cm_curve_plot_list_paad_niche"]])
saveRDS(res.cox.tcga.niche, snakemake@output[["coxph_summary_paad_niche"]])
saveRDS(clin.assoc.tcga.niche, snakemake@output[["clin_assoc_paad_niche"]])
saveRDS(clin.assoc.tcga.pval.niche, snakemake@output[["clin_assoc_paad_pval_niche"]])


# PanCuRx survival analysis
## load PanCuRx clinical data for plotting
pancurx.clin.for.plotting <- readRDS(snakemake@input[["clin_data_for_plotting_pancurx"]])

## GSVA for the signatures in TCGA
gsva.es <- readRDS(snakemake@input[["gsva_es_pancurx"]])

## niche enrichment scores - GSVA
niche.gsva.es <- gsva.es |>
        column_to_rownames("sample_id")

names(niche.gsva.es) <- str_split(names(niche.gsva.es), "\\.", simplify = T)[, 2]
names(niche.gsva.es) <- gsub(" Rep", "", names(niche.gsva.es))
# names(niche.gsva.es) <- gsub("-|, ", "_", names(niche.gsva.es))
names(niche.gsva.es) <- gsub(" ", "_", names(niche.gsva.es))

print("Checking if the columns in niche.factors are present in niche.gsva.es")
table(colnames(niche.factors) %in% names(niche.gsva.es))
print("Checking if the columns in niche.gsva.es are present in niche.factors")
table(names(niche.gsva.es) %in% colnames(niche.factors))

# niche.factors <- niche.factors[,names(niche.gsva.es)]
niche.gsva.es <- niche.gsva.es[, colnames(niche.factors)]
niche.gsva.es <- as.matrix(niche.gsva.es)

niche.gsva.es <- niche.gsva.es %*% t(niche.factors)
niche.gsva.es <- as.data.frame(niche.gsva.es) |>
	rownames_to_column("sample_id")
names(niche.gsva.es) <- gsub(" ", "_", names(niche.gsva.es))

## niche enrichment scores - marker expression mean values
head(pancurx.clin.for.plotting)
niche.sig.exprs.mean <- pancurx.clin.for.plotting |>
        select(1, ends_with("_mean")) |>
        column_to_rownames("sample_id")

names(niche.sig.exprs.mean) <- str_split(names(niche.sig.exprs.mean), "\\.", simplify = T)[, 2]
names(niche.sig.exprs.mean) <- gsub(" Rep", "", names(niche.sig.exprs.mean))
# names(niche.sig.exprs.mean) <- gsub("-|, ", "_", names(niche.sig.exprs.mean))
names(niche.sig.exprs.mean) <- gsub(" ", "_", names(niche.sig.exprs.mean))
names(niche.sig.exprs.mean) <- gsub("_mean", "", names(niche.sig.exprs.mean))

print("Checking if the columns in niche.factors are present in niche.sig.exprs.mean")
table(colnames(niche.factors) %in% names(niche.sig.exprs.mean))
print("Checking if the columns in niche.sig.exprs.mean are present in niche.factors")
table(names(niche.sig.exprs.mean) %in% colnames(niche.factors))

niche.sig.exprs.mean <- niche.sig.exprs.mean[, colnames(niche.factors)]
niche.sig.exprs.mean <- as.matrix(niche.sig.exprs.mean)

niche.sig.exprs.mean <- niche.sig.exprs.mean %*% t(niche.factors)
niche.sig.exprs.mean <- as.data.frame(niche.sig.exprs.mean) |>
        rownames_to_column("sample_id")
names(niche.sig.exprs.mean) <- gsub(" ", "_", names(niche.sig.exprs.mean))


## survival analysis
p_list <- list()
clin.assoc.list <- list()
cox.list <- list()
cox.gsva.list <- list()
res.cut.list <- list()
res.cut.gsva.list <- list()

for (niche in grep("Niche", colnames(niche.gsva.es), value = T)) {
        niche.gsva.es.subset <- niche.gsva.es |>
                select(c("sample_id", niche))
        names(niche.gsva.es.subset) <- c("sample_id", "gsva")

        niche.sig.exprs.mean.subset <- niche.sig.exprs.mean |>
                select(c("sample_id", niche))
        names(niche.sig.exprs.mean.subset) <- c("sample_id", "mean")

        # join niche gsva estimate and expression score to clinical data
        pancurx.clin.geneset <- inner_join(niche.gsva.es.subset, pancurx.clin.for.plotting, by = "sample_id")
        pancurx.clin.geneset <- inner_join(niche.sig.exprs.mean.subset, pancurx.clin.geneset, by = "sample_id")

        names(niche.gsva.es.subset) <- c("sample_id", paste0(niche, "_gsva"))
        names(niche.sig.exprs.mean.subset) <- c("sample_id", paste0(niche, "_mean"))
        pancurx.clin.for.plotting <- inner_join(niche.gsva.es.subset, pancurx.clin.for.plotting, by = "sample_id")
        pancurx.clin.for.plotting <- inner_join(niche.sig.exprs.mean.subset, pancurx.clin.for.plotting, by = "sample_id")

        # further remove samples with missing/not useful information
        pancurx.clin.geneset <- pancurx.clin.geneset[!is.na(pancurx.clin.geneset$OSFromSurgery), ]

        # model association with clinical variables
        clin.assoc.list[["ploidy"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + ploidy")),
                data = pancurx.clin.for.plotting
        ))$coefficients
        clin.assoc.list[["purity"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + purity")),
                data = pancurx.clin.for.plotting
        ))$coefficients
        # clin.assoc.list[["sex"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + Sex")), data = pancurx.clin.for.plotting))$coefficients
        clin.assoc.list[["stage"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + Clinical_Staging")),
                data = pancurx.clin.for.plotting |>
                        #filter(!grepl("^I$|III", Clinical_Staging))
                        #mutate(Clinical_Staging = gsub("A$|B$", "", Clinical_Staging))
                        filter(!grepl("^I$|^II$", Clinical_Staging))
        ))$coefficients
        # clin.assoc.list[["stage2"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + PathologicalStaging")),
        #         data = pancurx.clin.for.plotting |>
        #                 filter(!grepl("I - III", PathologicalStaging)) |>
        #                 mutate(PathologicalStaging = gsub("^IA$|^IB$", "I", PathologicalStaging)) |>
        #                 mutate(PathologicalStaging = gsub("^III$|^IV$", "III - IV", PathologicalStaging))
        # ))$coefficients
        clin.assoc.list[["grade"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + Histological_Grade")),
                data = pancurx.clin.for.plotting |>
                        filter(!grepl("N/A|Other|Unknown", Histological_Grade)) |>
                        mutate(Histological_Grade = plyr::mapvalues(Histological_Grade, 
                                                                    from = c("Well differentiated", "Moderately differentiated", "Poorly differentiated", "Undifferentiated"), 
                                                                    to = c("G1", "G2", "G3", "G4")))
        ))$coefficients
        clin.assoc.list[["adjuvant"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + AdjuvantTreatmentType")),
                data = pancurx.clin.for.plotting |>
                        filter(!grepl("Unknown", AdjuvantTreatmentType))
        ))$coefficients
        clin.assoc.list[["adjuvant_outcome"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + AdjuvantBest_Therapy_Response")),
                data = pancurx.clin.for.plotting |>
                        filter(!grepl("N/A|Unknown", AdjuvantBest_Therapy_Response))
        ))$coefficients
        # clin.assoc.list[["macrophage"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + CD68.PPC.Pos.Pix.Perc.Stroma")),
        #         data = pancurx.clin.for.plotting |>
        #                 filter(!grepl("Discrepancy", CD68.PPC.Pos.Pix.Perc.Stroma))
        # ))$coefficients
        clin.assoc.list[["subtype"]][[niche]] <- summary(lm(as.formula(paste0("`", niche, score_for_clin_assoc, "` ~ 0 + Moffitt.mod.y")),
                data = pancurx.clin.for.plotting
        ))$coefficients

        # fit CoxPH model
        cox1 <- coxph(Surv(OSFromSurgery, patient.vital_status) ~ mean, data = pancurx.clin.geneset)
        cox2 <- coxph(Surv(OSFromSurgery, patient.vital_status) ~ gsva, data = pancurx.clin.geneset)
        cox.list[[niche]] <- cox1
        cox.gsva.list[[niche]] <- cox2

        # find optimal cutpoints to categorize geneset loading
        res.cut <- surv_cutpoint(pancurx.clin.geneset,
                time = "OSFromSurgery", event = "patient.vital_status",
                variables = c("mean"), minprop = 0.5
        )
        res.cut.gsva <- surv_cutpoint(pancurx.clin.geneset,
                time = "OSFromSurgery", event = "patient.vital_status",
                variables = c("gsva"), minprop = 0.5
        )
        # summary(res.cut)
        res.cut.list[[niche]] <- res.cut
        res.cut.gsva.list[[niche]] <- res.cut.gsva

        # plot(res.cut.gsva, "niche_gsva", palette = "npg")

        res.cat <- surv_categorize(res.cut)
        res.cat.gsva <- surv_categorize(res.cut.gsva)
        # head(res.cat.gsva)

        # create survival curves
        sfit1 <- survfit(Surv(OSFromSurgery, patient.vital_status) ~ mean, data = res.cat)
        sfit2 <- survfit(Surv(OSFromSurgery, patient.vital_status) ~ gsva, data = res.cat.gsva)
        # summary(sfit, times=seq(0,365*5,365))

        p_list <- c(
                p_list,
                ggsurvplot_list(list(sfit1, sfit2),
                        data = list(res.cat, res.cat.gsva),
                        legend.title = list("exprs", "gsva"),
                        conf.int = TRUE, pval = TRUE, risk.table = TRUE, title = niche
                )
        )
}
rm(niche)

# Extract data from cox.list
univ_results <- lapply(cox.gsva.list, function(x) {
        x <- summary(x)
        p.value <- x$waldtest["pvalue"]
        wald.test <- x$waldtest["test"]
        beta <- x$coefficients[1, 1] # coefficient beta
        HR <- x$coefficients[1, 2] # exp(beta)
        HR.confint.lower <- x$conf.int[1, "lower .95"]
        HR.confint.upper <- x$conf.int[1, "upper .95"]
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
        rownames_to_column("niche")

res.cox.pancurx.niche <- res.cox


# Extract data from clin.assoc.list
holder.list <- lapply(clin.assoc.list, function(clin.assoc) {
        clin.assoc <- lapply(clin.assoc, function(x) x[, 3]) # 1 for estimate, 2 for std.error, 3 for t value, 4 for p value
        as.data.frame(do.call(rbind, clin.assoc)) %>% rownames_to_column("niche")
})
holder.list <- lapply(names(holder.list), function(clin.var) {
        df <- holder.list[[clin.var]]
        names(df) <- plyr::mapvalues(names(df), from = c("V1"), to = c(clin.var))
        df
})
names(holder.list) <- names(clin.assoc.list)

print(holder.list[["adjuvant_outcome"]])

clin.assoc.pancurx.niche <- Reduce(function(x, y) {
        left_join(x, y, by = "niche")
}, holder.list)
holder <- res.cox |>
        select(niche, HR) |>
        mutate(HR = log(HR))
names(holder) <- c("niche", "log(HR)")
clin.assoc.pancurx.niche <- left_join(holder, clin.assoc.pancurx.niche, by = "niche")

rm(holder.list, holder)

holder.list <- lapply(clin.assoc.list, function(clin.assoc) {
        clin.assoc <- lapply(clin.assoc, function(x) x[, 4])
        as.data.frame(do.call(rbind, clin.assoc)) %>% rownames_to_column("niche")
})
holder.list <- lapply(names(holder.list), function(clin.var) {
        df <- holder.list[[clin.var]]
        names(df) <- plyr::mapvalues(names(df), from = c("V1"), to = c(clin.var))
        df
})
names(holder.list) <- names(clin.assoc.list)

print(holder.list[["adjuvant_outcome"]])

clin.assoc.pancurx.pval.niche <- Reduce(function(x, y) {
        left_join(x, y, by = "niche")
}, holder.list)
holder <- res.cox |> select(niche, p.adj)
names(holder) <- c("niche", "log(HR)")
clin.assoc.pancurx.pval.niche <- left_join(holder, clin.assoc.pancurx.pval.niche, by = "niche")

rm(holder.list, holder)

# TO SAVE
saveRDS(niche.gsva.es, snakemake@output[["gsva_es_pancurx_niche"]])
saveRDS(niche.sig.exprs.mean, snakemake@output[["sig_exprs_mean_pancurx_niche"]])
saveRDS(pancurx.clin.for.plotting, snakemake@output[["clin_data_for_plotting_pancurx_niche"]])
saveRDS(p_list, snakemake@output[["cm_curve_plot_list_pancurx_niche"]])
saveRDS(res.cox.pancurx.niche, snakemake@output[["coxph_summary_pancurx_niche"]])
saveRDS(clin.assoc.pancurx.niche, snakemake@output[["clin_assoc_pancurx_niche"]])
saveRDS(clin.assoc.pancurx.pval.niche, snakemake@output[["clin_assoc_pancurx_pval_niche"]])



