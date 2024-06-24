# Load libraries
suppressPackageStartupMessages({
        library(magrittr)
        library(tidyverse)
        library(sjstats)
        library(tidyr)
        library(stringr)
        library(SingleCellExperiment)
        library(reticulate)
        #library(survminer)
        library(readxl)
})

# load Toronto expression data
toronto.mRNA <- readRDS(snakemake@input[["toronto_bulk_data"]])

# load TCGA clinical data from 2019 Cell paper
tcga.clin <- read_excel(snakemake@input[["tcga_clin_data"]],
        sheet = "TCGA-CDR",
        na = c("", "[Not Available]", "[Not Applicable]", "[Unknown]")
)
paad.clin <- tcga.clin |> filter(type == "PAAD")


paad.clin$ajcc_pathologic_tumor_stage <- plyr::mapvalues(paad.clin$ajcc_pathologic_tumor_stage,
        from = c("Stage III", "Stage IV"),
        to = c("Stage III/IV", "Stage III/IV")
)

## load TCGA PAAD Purist basal/classical calls
paad.purist.subtype <- read_csv(snakemake@input[["tcga_purist_subtype_calls"]])

paad.purist.subtype <- paad.purist.subtype |>
        mutate(bcr_patient_barcode = substr(ID, 1, 12))

paad.clin <- left_join(paad.clin, paad.purist.subtype, by = "bcr_patient_barcode")

## SAVE TIDY CLINICAL DATA for TCGA PAAD
write_tsv(paad.clin, snakemake@output[["tidy_paad_clin_data"]])


# load COMPASS survival data
compass.clin <- read_tsv(snakemake@input[["compass_clin_data"]])

compass.clin$patient_vital_status <- plyr::mapvalues(compass.clin$Alive,
        from = c("Alive", "Dead"),
        to = c(0, 1)
) |> as.numeric()
compass.clin$Moffitt_Sig <- compass.clin$`Moffitt Sig`
compass.clin$`Moffitt Sig` <- NULL
compass.clin$PCSI_ID <- compass.clin$`PCSI ID`
compass.clin$`PCSI ID` <- NULL

# load Toronto/PanCuRx clinical data
toronto.clin <- colData(toronto.mRNA) |> as.data.frame()
toronto.meta <- read_csv(snakemake@input[["toronto_meta_data"]])

pancurx.meta <- toronto.meta |> filter(Project == "ICGC_resected")

## clean up clinical data from Amy Zhang (shared by Tiak)
pancurx.meta <- pancurx.meta |>
	mutate(sample_id = tumour_DNAName) |>
 select(-tumour_collapsed_cram, -normal_reference_collapsed_cram, -celluloid_segments)
pancurx.meta$donor <- gsub("I", "I_", pancurx.meta$donor)
pancurx.meta$sex <- substring(pancurx.meta$sex, 1, 1)

## tidy up pancurx clinical data
pancurx.clin.new <- full_join(toronto.clin, pancurx.meta, by = "sample_id")
pancurx.clin.new <- pancurx.clin.new |> filter(grepl("Pa", sample_id))

pancurx.clin.new <- pancurx.clin.new |>
        mutate(PCSI_ID = ifelse(!is.na(PCSI_ID), PCSI_ID, donor)) |>
        mutate(OSFromSurgery = ifelse(!is.na(OSFromSurgery), OSFromSurgery, survival.days)) |>
        mutate(Deceased.x = ifelse(!is.na(Deceased.x), Deceased.x, Deceased.y)) |>
        mutate(Sex = ifelse(!is.na(Sex), Sex, sex)) |>
        mutate(ClinicalStaging = ifelse(!is.na(ClinicalStaging), ClinicalStaging, clinical_stage))

pancurx.clin.new <- pancurx.clin.new |>
        filter(grepl("_P$|_P_", sample_id)) |>
        filter(!is.na(OSFromSurgery)) |>
        filter(!grepl("Unknown|N/A", OSFromSurgery))

pancurx.clin.new <- pancurx.clin.new |>
        arrange(PCSI_ID) |>
        distinct(PCSI_ID, .keep_all = TRUE)

### more tidy up of pancurx clinical data
pancurx.clin.new$patient.vital_status <- str_split(pancurx.clin.new$Followup_DiseaseStatus, " ", simplify = T)[, 1]
print("PanCuRx vital status from Michael's data:")
table(pancurx.clin.new$patient.vital_status)
table(is.na(pancurx.clin.new$patient.vital_status))

print("PanCuRx vital status updated by Tiak's data:")
table(pancurx.clin.new$Deceased.x)
table(is.na(pancurx.clin.new$Deceased.x))

pancurx.clin.new$patient.vital_status <- pancurx.clin.new$Deceased.x
pancurx.clin.new <- pancurx.clin.new[pancurx.clin.new$patient.vital_status != "Unknown", ]
pancurx.clin.new$patient.vital_status <- plyr::mapvalues(pancurx.clin.new$patient.vital_status,
        from = c(FALSE, TRUE),
        to = c(0, 1)
) |> as.numeric()

pancurx.clin.new$OSFromClinicalDiagnosis <- as.numeric(pancurx.clin.new$OSFromClinicalDiagnosis)
pancurx.clin.new$OSFromSurgery <- as.numeric(pancurx.clin.new$OSFromSurgery)
pancurx.clin.new$DFSFromSurgery <- as.numeric(pancurx.clin.new$DFSFromSurgery)

pancurx.clin.new$Clinical_Staging <- plyr::mapvalues(pancurx.clin.new$ClinicalStaging,
        from = c("I - III", "IB or IIB", "IIA or IIB"),
        to = c("II", "II", "II")
)
pancurx.clin.new <- pancurx.clin.new |> mutate(Clinical_Staging = ifelse(Clinical_Staging == "Unknown", clinical_stage, Clinical_Staging))

pancurx.clin.new$Histological_Grade <- plyr::mapvalues(pancurx.clin.new$HistologicalGrade,
        from = c("Moderate to Poor", "Poorly to undifferentiated"),
        to = c("Poorly differentiated", "Undifferentiated")
)

# load IHC data for estimating macrophage content in pancurx and compass
pancurx.IHC <- read_csv(snakemake@input[["pancurx_ihc_data"]])

## add macrophage content to pancurx.clin.new and compass.clin
pancurx.clin.new <- left_join(pancurx.clin.new, pancurx.IHC |> select(PCSI_ID, CD68.PPC.Pos.Pix.Perc.Stroma), by = "PCSI_ID")
compass.clin <- left_join(compass.clin, pancurx.IHC |> select(PCSI_ID, CD68.PPC.Pos.Pix.Perc.Stroma), by = "PCSI_ID")

pancurx.clin.new <- pancurx.clin.new |>
        mutate(CD68_stroma_status = ifelse(CD68.PPC.Pos.Pix.Perc.Stroma > median(CD68.PPC.Pos.Pix.Perc.Stroma, na.rm = TRUE), "high", "low"))
compass.clin <- compass.clin |>
        mutate(CD68_stroma_status = ifelse(CD68.PPC.Pos.Pix.Perc.Stroma > median(CD68.PPC.Pos.Pix.Perc.Stroma, na.rm = TRUE), "high", "low"))

print("Number of samples with CD68 stroma status in PanCuRx:")
pancurx.clin.new$CD68_stroma_status |>
        is.na() |>
        table()
table(pancurx.clin.new$CD68_stroma_status)
print("Number of samples with CD68 stroma status in COMPASS:")
compass.clin$CD68_stroma_status |>
        is.na() |>
        table()
table(compass.clin$CD68_stroma_status)

# load basal/classical calls for estimating basal/classical subtype in pancurx and compass
pancurx.subtype <- read_tsv(snakemake@input[["toronto_subtype_calls"]])

names(pancurx.subtype) <- c("sample_id", "Moffitt.mod", "sscsubtype")
pancurx.subtype$PCSI_ID <- str_split(pancurx.subtype$sample_id, "_", simplify = T)[, c(1, 2)]
pancurx.subtype$PCSI_ID <- paste0(pancurx.subtype$PCSI_ID[, 1], "_", pancurx.subtype$PCSI_ID[, 2])

## add basal/classical subtype to pancurx.clin.new and compass.clin
pancurx.clin.new <- left_join(pancurx.clin.new, pancurx.subtype, by = "sample_id")
compass.clin <- left_join(compass.clin, pancurx.subtype, by = "PCSI_ID")

print("Number of samples with Moffitt basal/classical subtype in PanCuRx:")
pancurx.clin.new$Moffitt.mod.y |>
        is.na() |>
        table()
table(pancurx.clin.new$Moffitt.mod.y)
print("Number of samples with Moffitt basal/classical subtype in COMPASS:")
compass.clin$Moffitt.mod |>
        is.na() |>
        table()
table(compass.clin$Moffitt.mod)

print("Number of samples with Chan basal/classical subtype in PanCuRx:")
pancurx.clin.new$sscsubtype |>
        is.na() |>
        table()
table(pancurx.clin.new$sscsubtype)
print("Number of samples with Chan basal/classical subtype in COMPASS:")
compass.clin$sscsubtype |>
        is.na() |>
        table()
table(compass.clin$sscsubtype)

## tidy up basal/classical subtype calls in pancurx.clin.new and compass.clin
pancurx.clin.new <- pancurx.clin.new |>
        mutate(sscsubtype_coarse = gsub("A|B", "", sscsubtype))
pancurx.clin.new$sscsubtype_coarse <- gsub("asal", "Basal", pancurx.clin.new$sscsubtype_coarse)

compass.clin <- compass.clin |>
        mutate(sscsubtype_coarse = gsub("A|B", "", sscsubtype))
compass.clin$sscsubtype_coarse <- gsub("asal", "Basal", compass.clin$sscsubtype_coarse)

print("Number of samples with coarse Chan basal/classical subtype in PanCuRx:")
pancurx.clin.new$sscsubtype_coarse |>
        is.na() |>
        table()
table(pancurx.clin.new$sscsubtype_coarse)
print("Number of samples with coarse Chan basal/classical subtype in COMPASS:")
compass.clin$sscsubtype_coarse |>
        is.na() |>
        table()
table(compass.clin$sscsubtype_coarse)

## print basic information
print("Number of samples in PanCuRx:")
nrow(pancurx.clin.new)
print("Number of patients in PanCuRx:")
length(unique(pancurx.clin.new$PCSI_ID.x))

print("Number of samples in COMPASS:")
nrow(compass.clin)
print("Number of patients in COMPASS:")
length(unique(compass.clin$PCSI_ID))

# SAVE TIDY CLINICAL DATA for PanCuRx and COMPASS
write_tsv(pancurx.clin.new, snakemake@output[["tidy_pancurx_clin_data"]])
write_tsv(compass.clin, snakemake@output[["tidy_compass_clin_data"]])

# compare different versions of the clinical data of the Toronto cohort
print("Number of samples in both Michael's bulk data colData and Tiak's (from Amy Zhang, Jun2023) Toronto clinical data:")
intersect(toronto.clin$sample_id, toronto.meta$tumour_DNAName) |> length()
print("Number of samples in both Kieran's COMPASS_clean.tsv and Tiak's (from Amy Zhang, Jun2023) Toronto clinical data:")
intersect(compass.clin$sample_id, toronto.meta$tumour_DNAName) |> length()
print("Number of donors in both Kieran's COMPASS_clean.tsv and Tiak's (from Amy Zhang, Jun2023) Toronto clinical data:")
intersect(gsub("_", "", compass.clin$PCSI_ID), toronto.meta$donor) |> length()

#toronto.meta.test <- toronto.meta |> filter(donor %in% intersect(gsub("_", "", compass.clin$PCSI_ID), toronto.meta$donor))
toronto.meta.test <- toronto.meta |> filter(tumour_DNAName %in% intersect(compass.clin$sample_id, toronto.meta$tumour_DNAName))
toronto.meta.test$sample_id <- toronto.meta.test$tumour_DNAName

compass.clin.test <- compass.clin |> filter(sample_id %in% intersect(compass.clin$sample_id, toronto.meta$tumour_DNAName))
toronto.clin.test <- toronto.clin |> filter(sample_id %in% intersect(compass.clin$sample_id, toronto.meta$tumour_DNAName))

test.df <- full_join(toronto.clin.test |> select(sample_id, tissue, OSFromSurgery, DFSFromSurgery, Deceased, ClinicalStaging, PathologicalStaging, Stage),
        toronto.meta.test |> select(sample_id, survival.days, Deceased, clinical_stage),
        by = "sample_id",
        suffix = c("_oldmeta", "_toronto_Tiak")
)
test.df <- full_join(test.df,
        compass.clin.test |> select(sample_id, survival.days, Alive),
        by = "sample_id",
        suffix = c("_toronto_Tiak", "_compass_Kieran")
)
test.df <- test.df |>
        mutate(survival.days_agree = (survival.days_toronto_Tiak == survival.days_compass_Kieran)) |>
        mutate(Deceased_agree = (Deceased_oldmeta == Deceased_toronto_Tiak))

names(test.df) <- plyr::mapvalues(names(test.df),
        from = c(
                "tissue", "OSFromSurgery", "Alive",
                "ClinicalStaging", "PathologicalStaging", "clinical_stage"
        ),
        to = c(
                "tissue_oldmeta", "OSFromSurgery_oldmeta", "Alive_compass_Kieran",
                "ClinicalStaging_oldmeta", "PathologicalStaging_oldmeta", "clinical_stage_toronto_Tiak"
        )
)

test.df <- test.df |>
	select(sample_id, tissue_oldmeta,
	       OSFromSurgery_oldmeta, survival.days_compass_Kieran, survival.days_toronto_Tiak, survival.days_agree,
	       Deceased_oldmeta, Alive_compass_Kieran, Deceased_toronto_Tiak, Deceased_agree,
	       ClinicalStaging_oldmeta, PathologicalStaging_oldmeta, clinical_stage_toronto_Tiak)

write_tsv(test.df, snakemake@output[["toronto_compass_clin_data_comparison"]])
rm(toronto.meta.test, compass.clin.test, toronto.clin.test)

