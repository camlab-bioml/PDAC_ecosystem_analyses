suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(stats)
  library(gdata)
  library(magrittr)
  library(stringr)
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
})

# load stan data and other parameters
params_init_value_list <- readRDS(snakemake@input[["params_init_value_list"]])
fit.optim <- readRDS(snakemake@input[['fit_optim']])
stan.data <- readRDS(snakemake@input[["stan_data"]])
df.sig.mean <- readRDS(snakemake@input[["df_sig_mean"]])
sigs.encodings <- read_tsv(snakemake@input[["sigs_encodings"]])
samples.encodings <- read_tsv(snakemake@input[["samples_encodings"]])

celltypes <- snakemake@params[["celltypes"]]
number.of.niches <- snakemake@params[["number_of_niches"]] |> as.numeric()
nIter <- snakemake@params[["nIter"]] |> as.numeric()

# tidy up stan data data names
names(stan.data) <- gsub(" ", "_", names(stan.data))
names(stan.data) <- gsub("-", "_", names(stan.data))
names(stan.data) <- gsub(",", "", names(stan.data))

celltypes <- gsub(" ", "_", celltypes)
celltypes <- gsub("-", "_", celltypes)
celltypes <- gsub(",", "", celltypes)
print(celltypes)

#sigs.encodings$sig <- gsub(" ", "_", sigs.encodings$sig)
#sigs.encodings$sig <- gsub("-", "_", sigs.encodings$sig)
#sigs.encodings$sig <- gsub(",", "", sigs.encodings$sig)

# extract initial values
niche.factors.init.value <- params_init_value_list$niche_factors
niche.loadings.init.value <- params_init_value_list$niche_loadings

colnames(niche.factors.init.value) <- gsub("([0-9])$|([0-9][0-9])$", " \\1\\2", colnames(niche.factors.init.value))
rownames(niche.factors.init.value) <- paste("Niche ", seq_len(nrow(niche.factors.init.value)), sep = " ")
colnames(niche.loadings.init.value) <- paste("Niche ", seq_len(ncol(niche.loadings.init.value)), sep = " ")

saveRDS(niche.factors.init.value, file = snakemake@output[['microenvironment_niche_factors_init_value']])
saveRDS(niche.loadings.init.value, file = snakemake@output[['niche_factor_loadings_init_value']])

# extract parameters
holder <- fit.optim$mle(variables = "patient_specific_modelled_mu")
holder <- matrix(holder, nrow = stan.data$P)

holder1 <- fit.optim$mle(variables = "niche_factors")
holder1 <- matrix(holder1, nrow = stan.data$L)

holder2 <- fit.optim$mle(variables = "niche_loadings")
holder2 <- matrix(holder2, nrow = stan.data$P)

lapply(celltypes, function(ct) paste(ct, seq(stan.data[[paste0("K_", ct)]]), sep = " ")) %>% unlist()

colnames(holder) <- gsub(" Rep ", " ", gsub(" RepVal ", " ", sigs.encodings$sig))
rownames(holder) <- samples.encodings$sample
colnames(holder1) <- gsub(" Rep ", " ", gsub(" RepVal ", " ", sigs.encodings$sig))
rownames(holder1) <- paste("Niche ", seq_len(nrow(holder1)), sep = " ")
colnames(holder2) <- paste("Niche ", seq_len(ncol(holder2)), sep = " ")
rownames(holder2) <- samples.encodings$sample

saveRDS(holder, file = snakemake@output[['patient_specific_modelled_mu']])
saveRDS(holder1, file = snakemake@output[['microenvironment_niche_factors']])
saveRDS(holder2, file = snakemake@output[['niche_factor_loadings']])

holder3 <- lapply(celltypes, function(ct) {
  cov_i_mtx <- fit.optim$mle(variables = paste0("cov_i_", ct))
  cov_i_mtx <- matrix(cov_i_mtx, nrow = stan.data[[paste0("K_", ct)]])
  sig.names <- sigs.encodings |> filter(grepl(ct, sig)) |> pull(sig)
  colnames(cov_i_mtx) <- gsub(" Rep ", " ", gsub(" RepVal ", " ", sig.names))
  rownames(cov_i_mtx) <- gsub(" Rep ", " ", gsub(" RepVal ", " ", sig.names))
  cov_i_mtx
})
names(holder3) <- paste0("cov_i_", celltypes)

saveRDS(holder3, file = snakemake@output[['intrinsic_covariance_matrices']])