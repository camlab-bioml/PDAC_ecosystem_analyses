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

# set up cmdStanR
parallelly::availableCores()
options(mc.cores = 16)

check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
color_scheme_set("brightblue")
check_cmdstan_toolchain()

#install_cmdstan(cores = 2)

cmdstan_path()
cmdstan_version()

# load stan data and other parameters
stan.data <- readRDS(snakemake@input[["stan_data"]])
celltypes <- snakemake@params[["celltypes"]]
df.sig.mean <- readRDS(snakemake@input[["df_sig_mean"]])

number.of.niches <- snakemake@params[["number_of_niches"]] |> as.numeric()

nIter <- snakemake@params[["nIter"]] |> as.numeric()
nWarmup <- snakemake@params[["nWarmup"]]
nChains <- snakemake@params[["nChains"]]
nCores <- getOption("mc.cores", 1)
treeDepth <- snakemake@params[["treeDepth"]]

# tidy up stan data data names
names(stan.data) <- gsub(" ", "_", names(stan.data))
names(stan.data) <- gsub("-", "_", names(stan.data))
names(stan.data) <- gsub(",", "", names(stan.data))

celltypes <- gsub(" ", "_", celltypes)
celltypes <- gsub("-", "_", celltypes)
celltypes <- gsub(",", "", celltypes)

# load and compile model (cmdstanr)
model <- cmdstan_model(stan_file = snakemake@input[["stan_model"]])

model$exe_file()
model$print()

# run cmdstanr
## Optimization
### initial values for the variables
params.init.value.list <- lapply(celltypes, function(ct) {
  cov(stan.data[[paste0("y_", ct)]], method = "pearson")
})
names(params.init.value.list) <- paste0("cov_i_", celltypes)

# NMF is in conflict with DelayedArray because DelayedArray implemented a generic seed method as well
#saveRDS(df.sig.mean, "/Users/kieranlab/Desktop/Snakemake_pipelines/PDAC_TME/prototype/Stan/df_sig_mean_for_stan.rds")

nmf.niches <- NMF::nmf(x = as.matrix(df.sig.mean), rank = number.of.niches, method = "brunet", seed = NULL, model = NULL, nrun = 1, .options = 'p4tv')

#nmf.niches <- readRDS("/Users/kieranlab/Desktop/Snakemake_pipelines/PDAC_TME/prototype/Stan/nmf_niches_for_stan.rds")

params.init.value.list[["niche_loadings"]] <- nmf.niches@fit@W
params.init.value.list[["niche_factors"]] <- nmf.niches@fit@H

saveRDS(params.init.value.list, snakemake@output[['params_init_value_list']])

fit.optim <- model$optimize(
  data = stan.data,
  refresh = ceiling(nIter/100),
  init = list(
    params.init.value.list
  ),
  algorithm = "lbfgs",
  tol_param = 1e-7,
  iter = nIter
)

fit.optim$save_object(file = snakemake@output[['fit_optim']])
fit.optim$print()