suppressPackageStartupMessages({
        library(tidyverse)
        library(ComplexHeatmap)
})

# Load the data
## Discovery cohort
patient.profiles.df <- read_tsv(snakemake@input[["patient_profiles_collapsed"]])

### remove CD45 sorted samples, get meta data
patient.profiles.df <- filter(patient.profiles.df, !grepl("CD45m", sample))
patient.profiles.meta <- select(patient.profiles.df, cohort, sample)

write_tsv(patient.profiles.meta, snakemake@output[["patient_profiles_collapsed_meta"]])

### gather into tidy format and get a cell type for each signature
patient.profiles.df2 <- select(patient.profiles.df, -cohort) |>
  gather(signature, value, -sample) |>
  replace_na(list(value = 0)) |>
  mutate(cell_type = gsub("| Rep [0-9]$", "", signature)) |>
  mutate(cell_type = gsub(" Rep [0-9][0-9]$", "", cell_type))

table(patient.profiles.df2$cell_type)

### scale the signatures so they have a max value 1 and min value 0
patient.profiles.df2 <- group_by(patient.profiles.df2, cell_type, signature) |>
  mutate(minn = min(value, na.rm = TRUE),
         maxx = max(value, na.rm = TRUE)) |>
  mutate(scaled_value = (value - minn) / (maxx - minn)) |>
  #mutate(scaled_value = value) |>
  ungroup()

### un-tidy this back into matrix form
patient.profiles.df3 <- select(patient.profiles.df2, sample, signature, scaled_value) |> 
  spread(signature, scaled_value)

patient.profiles.mat <- select(patient.profiles.df3, -sample) |>
        as.matrix()
rownames(patient.profiles.mat) <- patient.profiles.df3$sample

patient.profiles.mat <- patient.profiles.mat[patient.profiles.meta$sample, ]

### save the cleaned and scaled matrix
patient.profiles.df.for.save <- patient.profiles.df3 |>
  mutate(cohort = plyr::mapvalues(patient.profiles.df3$sample,
				  from = patient.profiles.meta$sample,
				  to = patient.profiles.meta$cohort))

write_tsv(patient.profiles.df.for.save, snakemake@output[["patient_profiles_collapsed_cleaned_scaled"]])




## Validation cohort
patient.validation.profiles.df <- read_tsv(snakemake@input[["patient_profiles_collapsed_scored_val"]])

### get meta data
patient.validation.profiles.meta <- select(patient.validation.profiles.df, cohort, sample)

write_tsv(patient.validation.profiles.meta, snakemake@output[["patient_profiles_collapsed_scored_val_meta"]])

### normalize as per the discovery cohort
patient.validation.profiles.df2 <- select(patient.validation.profiles.df, -cohort) |>
        gather(signature, value, -sample) |>
        replace_na(list(value = 0)) |>
        mutate(cell_type = gsub("| RepVal [0-9]$", "", signature)) |>
        mutate(cell_type = gsub(" RepVal [0-9][0-9]$", "", cell_type)) |>
        group_by(cell_type, signature) |>
        mutate(
                minn = min(value, na.rm = TRUE),
                maxx = max(value, na.rm = TRUE)
        ) |>
        mutate(scaled_value = (value - minn) / (maxx - minn)) |>
        # mutate(scaled_value = value) |>
        ungroup()

table(patient.validation.profiles.df2$cell_type)

## un-tidy this back into matrix form
patient.validation.profiles.df3 <- select(patient.validation.profiles.df2, sample, signature, scaled_value) |>
        spread(signature, scaled_value)

patient.validation.profiles.mat <- select(patient.validation.profiles.df3, -sample) |>
        as.matrix()
rownames(patient.validation.profiles.mat) <- patient.validation.profiles.df3$sample

patient.validation.profiles.mat <- patient.validation.profiles.mat[patient.validation.profiles.meta$sample, ]

### save the cleaned and scaled matrix
patient.validation.profiles.df.for.save <- patient.validation.profiles.df3 |>
        mutate(cohort = plyr::mapvalues(patient.validation.profiles.df3$sample,
                from = patient.validation.profiles.meta$sample,
                to = patient.validation.profiles.meta$cohort
        ))

write_tsv(patient.validation.profiles.df.for.save, snakemake@output[["patient_profiles_collapsed_scored_val_cleaned_scaled"]])


# Are signature correlations consistent between cohorts?
colnames(patient.profiles.mat) <- gsub(" Rep", "", colnames(patient.profiles.mat))
colnames(patient.validation.profiles.mat) <- gsub(" RepVal", "", colnames(patient.validation.profiles.mat))

cc.disc <- cor(patient.profiles.mat, method = "spearman")
cc.val <- cor(patient.validation.profiles.mat, method = "spearman")
cc <- cc.disc

dfc.disc <- as.data.frame(cc.disc) |>
        rownames_to_column("signature_1") |>
        gather(signature_2, correlation_discovery, -signature_1)

dfc.val <- as.data.frame(cc.val) |>
        rownames_to_column("signature_1") |>
        gather(signature_2, correlation_validation, -signature_1)

head(dfc.disc)
head(dfc.val)

## join these together
dfc <- inner_join(dfc.disc, dfc.val) |>
        as_tibble() |>
        filter(signature_1 != signature_2)

dfc

## are within-cell type correlations more correlated than between?
dfc <- dfc |>
        mutate(cell_type_1 = gsub("| [0-9]$", "", signature_1)) |>
        mutate(cell_type_2 = gsub("| [0-9]$", "", signature_2)) |>
        mutate(cell_type_1 = gsub("| [0-9][0-9]$", "", cell_type_1)) |>
        mutate(cell_type_2 = gsub("| [0-9][0-9]$", "", cell_type_2)) |>
        mutate(same_cell_type = cell_type_1 == cell_type_2)

## let's not double-count correlations
dfc$cell_type_str <- apply(dfc, 1, function(x) {
        x <- x[c("signature_1", "signature_2")]
        x <- sort(x)
        paste(x, collapse = "_")
})

dfc <- dfc[!duplicated(dfc$cell_type_str), ]

head(dfc)

### save to patient profile correlation data frame
write_tsv(dfc, snakemake@output[["patient_profiles_correlation_data_frame"]])

with(dfc, cor.test(correlation_discovery, correlation_validation))

with(filter(dfc, same_cell_type), cor.test(correlation_discovery, correlation_validation))

with(filter(dfc, !same_cell_type), cor.test(correlation_discovery, correlation_validation))
