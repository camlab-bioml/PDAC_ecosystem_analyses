functions {
  
}

data {
  int<lower = 1> C;    // number of cell types
  
  int<lower = 1> L; // number of 'niches'

  int<lower = 1> N_mono;   // number of monocyte cells 
  int<lower = 1> N_ste;   // number of activated stellate cells 

  int<lower = 1> P_mono;   // number of patients having monocyte signatures
  int<lower = 1> P_ste;   // number of patients having activated stellate signatures

  int<lower = 1> K_mono;   // number of monocyte signatures
  int<lower = 1> K_ste;    // number of activated stellate signatures

  matrix[N_mono, K_mono] y_mono;     // monocyte signature loading matrix
  matrix[N_ste, K_ste] y_ste;     // activated stellate signature loading matrix

  int<lower = 1, upper = P_mono> x_mono[N_mono];    // patient id for monocytes single-cells 
  int<lower = 1, upper = P_ste> x_ste[N_ste];    // patient id for activated stellate single-cells 
}

transformed data {
  int<lower = 1> K;
  int<lower = 1> P;
  K = K_mono + K_ste;  


  P = P_mono;
}

parameters {
  
  // patient specific monocyte means
  matrix<lower = 0>[P_mono, K_mono] mu_mono;
  // monocyte intrinic variance
  cov_matrix[K_mono] cov_i_mono;
  
  matrix<lower = 0>[P_ste, K_ste] mu_ste;
  cov_matrix[K_ste] cov_i_ste;
  
  // vector<lower=0>[K_mono] sigma_mono; // standard deviation of each signature

  matrix<lower=0> [L,K]  niche_factors;
  matrix<lower=0> [P,L]  niche_loadings;

}

transformed parameters {
  matrix<lower = 0>[P, K] theta;
  matrix[P,K] patient_specific_modelled_mu;
  
  theta = append_col(mu_mono, mu_ste);
  patient_specific_modelled_mu = niche_loadings * niche_factors;


}

model {
  
  // priors on intrinsic covariance matrices
  target += inv_wishart_lpdf(cov_i_mono | K_mono, 0.01 + diag_matrix(rep_vector(1.0,K_mono)));
  target += inv_wishart_lpdf(cov_i_ste | K_ste, 0.01 + diag_matrix(rep_vector(1.0,K_ste)));
  
  // priors on niche factors
  for(l in 1:L) {
    target += normal_lpdf(niche_factors[l] | 0, 10);
  }
  
  for(p in 1:P) {
    target += normal_lpdf(niche_loadings[p] | 0, 1);
    target += normal_lpdf(theta[p] | patient_specific_modelled_mu[p], 1);
  }
  
  for(n in 1:N_mono) {
    int patient_n = x_mono[n];
    target += multi_normal_lpdf(y_mono[n] | mu_mono[patient_n], cov_i_mono);
  }
  
  for(n in 1:N_ste) {
    int patient_n = x_ste[n];
    target += multi_normal_lpdf(y_ste[n] | mu_ste[patient_n], cov_i_ste);
  }
  
//   
//   //prior 
//   for (p in P_list) {
//   	// target += normal_lpdf(mu_mono[p] | 0, 10);
//   	// target += normal_lpdf(mu_act_ste[p] | 0, 10);
//   }
// 
// //   target += inv_wishart_lpdf(cov_e_mono | K_mono, diag_matrix(rep_vector(1.0, K_mono)));
// //   target += inv_wishart_lpdf(cov_e_act_ste | K_act_ste, diag_matrix(rep_vector(1.0, K_act_ste)));
// //   target += inv_wishart_lpdf(cov_i_mono | K_mono, diag_matrix(rep_vector(1.0, K_mono)));
// //   target += inv_wishart_lpdf(cov_i_act_ste | K_act_ste, diag_matrix(rep_vector(1.0, K_act_ste)));
// 
//   sigma_mono ~ cauchy(0, 5); // prior for sigma
//   Lcorr_mono ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
//   sigma_act_ste ~ cauchy(0, 5); // prior for sigma
//   Lcorr_act_ste ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
// 
//   target += normal_lpdf(gamma | 0, 1);
//   target += inv_wishart_lpdf(cov_e | K, diag_matrix(rep_vector(1.0, K)));
// 
// 
//   //likelihood
// //   for (p in 1:P_mono) {
// // 	target += multi_normal_lpdf(mu_mono[p] | eta_mono, cov_e_mono); 
// //   }
//   for (n in 1:N_mono) {
// 	target += multi_normal_lpdf(y_mono[n] | mu_mono[x_mono[n]], cov_i_mono);
//   }
//   
// //   for (p in 1:P_act_ste) {
// // 	target += multi_normal_lpdf(mu_act_ste[p] | eta_act_ste, cov_e_act_ste); 
// //   }
//   for (n in 1:N_act_ste) {
// 	target += multi_normal_lpdf(y_act_ste[n] | mu_act_ste[x_act_ste[n]], cov_i_act_ste);
//   }
// 
//   //likelihood
//   for (p in P_list) {
//      target += multi_normal_lpdf(theta[p] | gamma, cov_e);
//   }
  
}

generated quantities {
//   matrix[N_mono, K_mono] Y_mono;     // monocyte signature loading matrix
//   matrix[N_act_ste, K_act_ste] Y_act_ste;     // activated stellate signature loading matrix
// 
//   for (n in 1:N_mono) {
// 	Y_mono[n] = to_row_vector(multi_normal_rng(mu_mono[x_mono[n]], cov_i_mono));
//   }
//   for (n in 1:N_act_ste) {
// 	Y_act_ste[n] = to_row_vector(multi_normal_rng(mu_act_ste[x_act_ste[n]], cov_i_act_ste));
//   }

} // The posterior predictive distribution