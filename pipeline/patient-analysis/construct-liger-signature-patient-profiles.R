suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(sjstats)
  library(tidyr)
  library(stringr)
})

profile.flavor <- snakemake@wildcards[['profile']]
condition <- snakemake@wildcards[['condition']]
compartment <- snakemake@wildcards[['compartment']]

celltypes <- str_split(snakemake@input[['celltype_sig_profiles']], "signature-analysis/", simplify = T)[,2]
celltypes <- str_split(celltypes, "/signature-loading-profiles", simplify = T)[,1]

# load signature loading profiles for each cell type
sigprofile.dflist <- lapply(snakemake@input[['celltype_sig_profiles']], read_tsv)
names(sigprofile.dflist) <- celltypes

for (ct in celltypes) {
  colnames(sigprofile.dflist[[ct]]) <- gsub(condition, ct, colnames(sigprofile.dflist[[ct]]))
}
rm(ct)

# combine signature loading profiles from all cell types
sigprofile.df <- Reduce(function(df1, df2) full_join(df1, df2, by = c("cohort", "sample")), sigprofile.dflist)

# long form patient profiles
sigprofile.df.long <- sigprofile.df %>%
  pivot_longer(cols = !contains("cohort") & !contains("sample"), names_to = "signature", values_to = profile.flavor) %>%
  mutate(celltype = str_split(signature, pattern = " ", simplify = T)[,1])

# get simple correlation across signatures
sigprofile.mtx <- sigprofile.df %>% 
  select(!c(cohort, sample)) %>%
  as.matrix()

print("Summarizing all signatures: ")
summary(sigprofile.mtx %>% c())

sigprofile.corr <- cor(sigprofile.mtx, use = "pairwise.complete.obs")

# save the results
write_tsv(sigprofile.df, file = snakemake@output[['patient_profiles']])
write_tsv(sigprofile.df.long, file = snakemake@output[['patient_profiles_long']])
write_tsv(as.data.frame(sigprofile.corr), file = snakemake@output[['patient_profiles_corr']])



