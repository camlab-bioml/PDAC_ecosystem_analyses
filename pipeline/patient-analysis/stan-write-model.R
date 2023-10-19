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

install_cmdstan(cores = 2)

cmdstan_path()
cmdstan_version()

# load stan data and other parameters
stan.data <- readRDS(snakemake@input[["stan_data"]])
celltypes <- snakemake@params[["celltypes"]]
lambda <- snakemake@params[["lambda_for_orthogonality"]]

# tidy up stan data data names
names(stan.data) <- gsub(" ", "_", names(stan.data))
names(stan.data) <- gsub("-", "_", names(stan.data))
names(stan.data) <- gsub(",", "", names(stan.data))

# cmdstanr
## write model
celltypes <- gsub(" ", "_", celltypes)
celltypes <- gsub("-", "_", celltypes)
celltypes <- gsub(",", "", celltypes)

stan_code <- paste0(
"
functions {
  vector get_upper_triangular(matrix m) {
    int K = num_elements(m);
    print('K = ', K);
    int L = size(m);
    print('L = ', L);
    vector[K-(L*(L+1)%/%2)] v;
    int idx = 1;
    for (i in 1:rows(m)) {
      for (j in 1:cols(m)) {
        if (j > i) {
          v[idx] = m[i, j];
          idx = idx + 1;
        }
      }
    }
    return v;
  }
}

data {
  int<lower = 1> C;    // number of cell types
  
  int<lower = 1> L; // number of 'niches'
  
  int<lower = 1> P; // number of patients
  
", 
(paste0("
  int<lower = 1> j_", celltypes, ";   // starting column index for ", celltypes, " cells 
  int<lower = 1> N_", celltypes, ";   // number of ", celltypes, " cells 
  int<lower = 1> P_", celltypes, ";   // number of ", celltypes, " patients
  int<lower = 1> K_", celltypes, ";   // number of ", celltypes, " signatures
  matrix[N_", celltypes, ", K_", celltypes, "] y_", celltypes, ";     // ", celltypes, " signature loading matrix
  array[N_", celltypes, "] int<lower = 1, upper = P_", celltypes, "> x_", celltypes, ";    // patient id for ", celltypes, " single-cells
") %>% str_flatten()), 
"
}

transformed data {
  int<lower = 1> K;
  // int<lower = 1> P;
  
  K = 0", 
(paste0(" + K_", celltypes
) %>% str_flatten()), ";",
(paste0("
  // P = P_", celltypes[1], ";
") %>% str_flatten()),
"}

parameters {
  // intrinsic variance
", 
(paste0("  cov_matrix[K_", celltypes, "] cov_i_", celltypes, ";
") %>% str_flatten()), 
"
  matrix<lower = 0> [P, K] theta;
  matrix<lower = 0> [P, L] niche_loadings;
  matrix<lower = 0> [L, K] niche_factors;
}

transformed parameters {
  matrix[L, L] niche_factors_prod;
  niche_factors_prod = tcrossprod(niche_factors); // tcrossprod(x) = xx' = xx^T

  vector<lower = 0>[L*L - (L*(L+1)%/%2)] niche_factors_prod_upper_triangular;
  niche_factors_prod_upper_triangular = get_upper_triangular(niche_factors_prod);
  
  matrix[P,K] patient_specific_modelled_mu;
  patient_specific_modelled_mu = niche_loadings * niche_factors;
  
  // patient specific means
", 
(paste0("
  matrix<lower = 0>[P, K_", celltypes, "] mu_", celltypes, ";
  mu_", celltypes, " = block(theta, 1, j_", celltypes, ", P, K_", celltypes, ");
") %>% str_flatten()), 
"  
}

model {
  
  // priors on intrinsic covariance matrices
  ", 
(paste0("target += inv_wishart_lpdf(cov_i_", celltypes, " | K_", celltypes, ", 0.01 + diag_matrix(rep_vector(1.0, K_", celltypes, ")));
  ") %>% str_flatten()), 
"
  // priors on niche factors
  for(l in 1:L) {
    target += normal_lpdf(niche_factors[l] | 0, 10);
  }
  
  for(p in 1:P) {
    target += normal_lpdf(niche_loadings[p] | 0, 1);
    target += normal_lpdf(theta[p] | patient_specific_modelled_mu[p], 1);
  }
", 
(paste0("
  for(n in 1:N_", celltypes, ") {
    int patient_n = x_", celltypes, "[n];
    target += multi_normal_lpdf(y_", celltypes, "[n] | mu_", celltypes, "[patient_n], cov_i_", celltypes, ");
  }
") %>% str_flatten()), 
" 
  
  // orthogonality constraint
  target += -(", lambda, " * mean(niche_factors_prod_upper_triangular));

}

generated quantities {

",
(paste0("
  // for(n in 1:N_", celltypes, ") {
  //   int patient_n = x_", celltypes, "[n];
  //  y_", celltypes, "[n] = multi_normal_rng(mu_", celltypes, "[patient_n], cov_i_", celltypes, ");
  //}
") %>% str_flatten()),
"
} // The posterior predictive distribution
"
) 

write_stan_file(code = stan_code, dir = snakemake@params[["stan_model_dir"]], basename = snakemake@params[["stan_model_basename"]])
list.files(path = snakemake@params[["stan_model_dir"]])

## compile model (cmdstanr)
file <- paste0(snakemake@params[['stan_model_dir']], snakemake@params[['stan_model_basename']], ".stan")
model <- cmdstan_model(stan_file = file)

model$exe_file()
model$print()