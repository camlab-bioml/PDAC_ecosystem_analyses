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
fit.optim <- readRDS(snakemake@input[['fit_optim']])
stan.data <- readRDS(snakemake@input[["stan_data"]])
df.sig.mean <- readRDS(snakemake@input[["df_sig_mean"]])

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

# extract parameters
holder <- fit.optim$mle(variables = "patient_specific_modelled_mu")
holder <- matrix(holder, nrow = stan.data$P)

holder1 <- fit.optim$mle(variables = "niche_factors")
holder1 <- matrix(holder1, nrow = stan.data$L)

holder2 <- fit.optim$mle(variables = "niche_loadings")
holder2 <- matrix(holder2, nrow = stan.data$P)

lapply(celltypes, function(ct) paste(ct, seq(stan.data[[paste0("K_", ct)]]), sep = " ")) %>% unlist()

colnames(holder) <- lapply(celltypes, function(ct) paste(ct, seq(stan.data[[paste0("K_", ct)]]), sep = " ")) %>% unlist()
#rownames(holder) <- paste("Patient", seq_len(nrow(holder)), sep = " ")
colnames(holder1) <- lapply(celltypes, function(ct) paste(ct, seq(stan.data[[paste0("K_", ct)]]), sep = " ")) %>% unlist()
rownames(holder1) <- paste("Niche ", seq_len(nrow(holder1)), sep = " ")
colnames(holder2) <- paste("Niche ", seq_len(ncol(holder2)), sep = " ")
#rownames(holder2) <- paste("Patient", seq_len(nrow(holder2)), sep = " ")

saveRDS(holder, file = snakemake@output[['patient_specific_modelled_mu']])
saveRDS(holder1, file = snakemake@output[['microenvironment_niche_factors']])
saveRDS(holder2, file = snakemake@output[['niche_factor_loadings']])

holder3 <- lapply(celltypes, function(ct) {
  cov_i_mtx <- fit.optim$mle(variables = paste0("cov_i_", ct))
  cov_i_mtx <- matrix(cov_i_mtx, nrow = stan.data[[paste0("K_", ct)]])
  colnames(cov_i_mtx) <- paste(ct, seq(stan.data[[paste0("K_", ct)]]), sep = " ")
  rownames(cov_i_mtx) <- paste(ct, seq(stan.data[[paste0("K_", ct)]]), sep = " ")
  cov_i_mtx
})
names(holder3) <- paste0("cov_i_", celltypes)

saveRDS(holder3, file = snakemake@output[['intrinsic_covariance_matrices']])