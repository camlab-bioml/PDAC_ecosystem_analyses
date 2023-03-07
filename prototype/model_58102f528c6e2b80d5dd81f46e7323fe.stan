
functions {
  
}

data {
  int<lower = 1> C;    // number of cell types
  
  int<lower = 1> L; // number of 'niches'
  
  int<lower = 1> P; // number of patients

  int<lower = 1> j_Mono;   // starting column index for Mono cells 
  int<lower = 1> N_Mono;   // number of Mono cells 
  int<lower = 1> P_Mono;   // number of Mono patients
  int<lower = 1> K_Mono;   // number of Mono signatures
  matrix[N_Mono, K_Mono] y_Mono;     // Mono signature loading matrix
  array[N_Mono] int<lower = 1, upper = P_Mono> x_Mono;    // patient id for Mono single-cells

  int<lower = 1> j_stellate;   // starting column index for stellate cells 
  int<lower = 1> N_stellate;   // number of stellate cells 
  int<lower = 1> P_stellate;   // number of stellate patients
  int<lower = 1> K_stellate;   // number of stellate signatures
  matrix[N_stellate, K_stellate] y_stellate;     // stellate signature loading matrix
  array[N_stellate] int<lower = 1, upper = P_stellate> x_stellate;    // patient id for stellate single-cells

  int<lower = 1> j_ductal;   // starting column index for ductal cells 
  int<lower = 1> N_ductal;   // number of ductal cells 
  int<lower = 1> P_ductal;   // number of ductal patients
  int<lower = 1> K_ductal;   // number of ductal signatures
  matrix[N_ductal, K_ductal] y_ductal;     // ductal signature loading matrix
  array[N_ductal] int<lower = 1, upper = P_ductal> x_ductal;    // patient id for ductal single-cells

  int<lower = 1> j_CD8;   // starting column index for CD8 cells 
  int<lower = 1> N_CD8;   // number of CD8 cells 
  int<lower = 1> P_CD8;   // number of CD8 patients
  int<lower = 1> K_CD8;   // number of CD8 signatures
  matrix[N_CD8, K_CD8] y_CD8;     // CD8 signature loading matrix
  array[N_CD8] int<lower = 1, upper = P_CD8> x_CD8;    // patient id for CD8 single-cells

  int<lower = 1> j_B;   // starting column index for B cells 
  int<lower = 1> N_B;   // number of B cells 
  int<lower = 1> P_B;   // number of B patients
  int<lower = 1> K_B;   // number of B signatures
  matrix[N_B, K_B] y_B;     // B signature loading matrix
  array[N_B] int<lower = 1, upper = P_B> x_B;    // patient id for B single-cells

  int<lower = 1> j_CD4;   // starting column index for CD4 cells 
  int<lower = 1> N_CD4;   // number of CD4 cells 
  int<lower = 1> P_CD4;   // number of CD4 patients
  int<lower = 1> K_CD4;   // number of CD4 signatures
  matrix[N_CD4, K_CD4] y_CD4;     // CD4 signature loading matrix
  array[N_CD4] int<lower = 1, upper = P_CD4> x_CD4;    // patient id for CD4 single-cells

  int<lower = 1> j_DC;   // starting column index for DC cells 
  int<lower = 1> N_DC;   // number of DC cells 
  int<lower = 1> P_DC;   // number of DC patients
  int<lower = 1> K_DC;   // number of DC signatures
  matrix[N_DC, K_DC] y_DC;     // DC signature loading matrix
  array[N_DC] int<lower = 1, upper = P_DC> x_DC;    // patient id for DC single-cells

  int<lower = 1> j_endothelial;   // starting column index for endothelial cells 
  int<lower = 1> N_endothelial;   // number of endothelial cells 
  int<lower = 1> P_endothelial;   // number of endothelial patients
  int<lower = 1> K_endothelial;   // number of endothelial signatures
  matrix[N_endothelial, K_endothelial] y_endothelial;     // endothelial signature loading matrix
  array[N_endothelial] int<lower = 1, upper = P_endothelial> x_endothelial;    // patient id for endothelial single-cells

}

transformed data {
  int<lower = 1> K;
  #int<lower = 1> P;
  
  K = 0 + K_Mono + K_stellate + K_ductal + K_CD8 + K_B + K_CD4 + K_DC + K_endothelial;
  #P = P_Mono;
}

parameters {
  // intrinsic variance
  cov_matrix[K_Mono] cov_i_Mono;
  cov_matrix[K_stellate] cov_i_stellate;
  cov_matrix[K_ductal] cov_i_ductal;
  cov_matrix[K_CD8] cov_i_CD8;
  cov_matrix[K_B] cov_i_B;
  cov_matrix[K_CD4] cov_i_CD4;
  cov_matrix[K_DC] cov_i_DC;
  cov_matrix[K_endothelial] cov_i_endothelial;

  matrix<lower = 0> [P, K] theta;
  matrix<lower = 0> [P, L] niche_loadings;
  matrix<lower = 0> [L, K] niche_factors;
}

transformed parameters {
  matrix[P,K] patient_specific_modelled_mu;
  patient_specific_modelled_mu = niche_loadings * niche_factors;
  
  // patient specific means

  matrix<lower = 0>[P, K_Mono] mu_Mono;
  mu_Mono = block(theta, 1, j_Mono, P, K_Mono);

  matrix<lower = 0>[P, K_stellate] mu_stellate;
  mu_stellate = block(theta, 1, j_stellate, P, K_stellate);

  matrix<lower = 0>[P, K_ductal] mu_ductal;
  mu_ductal = block(theta, 1, j_ductal, P, K_ductal);

  matrix<lower = 0>[P, K_CD8] mu_CD8;
  mu_CD8 = block(theta, 1, j_CD8, P, K_CD8);

  matrix<lower = 0>[P, K_B] mu_B;
  mu_B = block(theta, 1, j_B, P, K_B);

  matrix<lower = 0>[P, K_CD4] mu_CD4;
  mu_CD4 = block(theta, 1, j_CD4, P, K_CD4);

  matrix<lower = 0>[P, K_DC] mu_DC;
  mu_DC = block(theta, 1, j_DC, P, K_DC);

  matrix<lower = 0>[P, K_endothelial] mu_endothelial;
  mu_endothelial = block(theta, 1, j_endothelial, P, K_endothelial);
  
}

model {
  
  // priors on intrinsic covariance matrices
  target += inv_wishart_lpdf(cov_i_Mono | K_Mono, 0.01 + diag_matrix(rep_vector(1.0, K_Mono)));
  target += inv_wishart_lpdf(cov_i_stellate | K_stellate, 0.01 + diag_matrix(rep_vector(1.0, K_stellate)));
  target += inv_wishart_lpdf(cov_i_ductal | K_ductal, 0.01 + diag_matrix(rep_vector(1.0, K_ductal)));
  target += inv_wishart_lpdf(cov_i_CD8 | K_CD8, 0.01 + diag_matrix(rep_vector(1.0, K_CD8)));
  target += inv_wishart_lpdf(cov_i_B | K_B, 0.01 + diag_matrix(rep_vector(1.0, K_B)));
  target += inv_wishart_lpdf(cov_i_CD4 | K_CD4, 0.01 + diag_matrix(rep_vector(1.0, K_CD4)));
  target += inv_wishart_lpdf(cov_i_DC | K_DC, 0.01 + diag_matrix(rep_vector(1.0, K_DC)));
  target += inv_wishart_lpdf(cov_i_endothelial | K_endothelial, 0.01 + diag_matrix(rep_vector(1.0, K_endothelial)));
  
  // priors on niche factors
  for(l in 1:L) {
    target += normal_lpdf(niche_factors[l] | 0, 10);
  }
  
  for(p in 1:P) {
    target += normal_lpdf(niche_loadings[p] | 0, 1);
    target += normal_lpdf(theta[p] | patient_specific_modelled_mu[p], 1);
  }

  for(n in 1:N_Mono) {
    int patient_n = x_Mono[n];
    target += multi_normal_lpdf(y_Mono[n] | mu_Mono[patient_n], cov_i_Mono);
  }

  for(n in 1:N_stellate) {
    int patient_n = x_stellate[n];
    target += multi_normal_lpdf(y_stellate[n] | mu_stellate[patient_n], cov_i_stellate);
  }

  for(n in 1:N_ductal) {
    int patient_n = x_ductal[n];
    target += multi_normal_lpdf(y_ductal[n] | mu_ductal[patient_n], cov_i_ductal);
  }

  for(n in 1:N_CD8) {
    int patient_n = x_CD8[n];
    target += multi_normal_lpdf(y_CD8[n] | mu_CD8[patient_n], cov_i_CD8);
  }

  for(n in 1:N_B) {
    int patient_n = x_B[n];
    target += multi_normal_lpdf(y_B[n] | mu_B[patient_n], cov_i_B);
  }

  for(n in 1:N_CD4) {
    int patient_n = x_CD4[n];
    target += multi_normal_lpdf(y_CD4[n] | mu_CD4[patient_n], cov_i_CD4);
  }

  for(n in 1:N_DC) {
    int patient_n = x_DC[n];
    target += multi_normal_lpdf(y_DC[n] | mu_DC[patient_n], cov_i_DC);
  }

  for(n in 1:N_endothelial) {
    int patient_n = x_endothelial[n];
    target += multi_normal_lpdf(y_endothelial[n] | mu_endothelial[patient_n], cov_i_endothelial);
  }
 
}

generated quantities {

} // The posterior predictive distribution

