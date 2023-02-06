functions {
  
}

data {
  int<lower = 1> C;    // number of cell types

  int<lower = 1> N_mono;   // number of monocyte cells 
  int<lower = 1> N_act_ste;   // number of activated stellate cells 

  int<lower = 1> P_mono;   // number of patients having monocyte signatures
  int<lower = 1> P_act_ste;   // number of patients having activated stellate signatures

  array[P_mono] int<lower = 1, upper = P_mono> P_list;

  int<lower = 1> K_mono;   // number of monocyte signatures
  int<lower = 1> K_act_ste;    // number of activated stellate signatures

  matrix[N_mono, K_mono] y_mono;     // monocyte signature loading matrix
  matrix[N_act_ste, K_act_ste] y_act_ste;     // activated stellate signature loading matrix

  array[N_mono] int<lower = 1, upper = P_mono> x_mono;    // patient id for monocytes single-cells 
  array[N_act_ste] int<lower = 1, upper = P_act_ste> x_act_ste;    // patient id for activated stellate single-cells 
}

transformed data {
  int<lower = 1> K;
  K = K_mono + K_act_ste;  

  int<lower = 1> P;
  P = P_mono;
}

parameters {
//   vector<lower = 0>[K_mono] eta_mono;	
  matrix<lower = 0>[P_mono, K_mono] mu_mono;
//   cov_matrix[K_mono] cov_e_mono;
//   cov_matrix[K_mono] cov_i_mono;
  cholesky_factor_corr[K_mono] Lcorr_mono; // cholesky factor (L_u matrix for R)
  vector<lower=0>[K_mono] sigma_mono; // standard deviation of each signature

//   vector<lower = 0>[K_act_ste] eta_act_ste;
  matrix<lower = 0>[P_act_ste, K_act_ste] mu_act_ste;
//   cov_matrix[K_act_ste] cov_e_act_ste;
//   cov_matrix[K_act_ste] cov_i_act_ste;
  cholesky_factor_corr[K_act_ste] Lcorr_act_ste; // cholesky factor (L_u matrix for R)
  vector<lower=0>[K_act_ste] sigma_act_ste; // standard deviation of each signature

  vector<lower = 0>[K] gamma;
  cov_matrix[K] cov_e;
}

transformed parameters {
  matrix<lower = 0>[P, K] theta;
  theta = append_col(mu_mono, mu_act_ste);

  corr_matrix[K_mono] R_mono; // correlation matrix
  cov_matrix[K_mono] cov_i_mono; // VCV matrix
  R_mono = multiply_lower_tri_self_transpose(Lcorr_mono); // R = Lcorr * Lcorr'
  cov_i_mono = quad_form_diag(R_mono, sigma_mono); // quad_form_diag: diag_matrix(sigma) * R * diag_matrix(sigma)

  corr_matrix[K_act_ste] R_act_ste; // correlation matrix
  cov_matrix[K_act_ste] cov_i_act_ste; // VCV matrix
  R_act_ste = multiply_lower_tri_self_transpose(Lcorr_act_ste); // R = Lcorr * Lcorr'
  cov_i_act_ste = quad_form_diag(R_act_ste, sigma_act_ste); // quad_form_diag: diag_matrix(sigma) * R * diag_matrix(sigma)
}

model {
  
  //prior 
  for (p in P_list) {
	target += normal_lpdf(mu_mono[p] | 0, 1);
	target += normal_lpdf(mu_act_ste[p] | 0, 1);
  }

//   target += inv_wishart_lpdf(cov_e_mono | K_mono, diag_matrix(rep_vector(1.0, K_mono)));
//   target += inv_wishart_lpdf(cov_e_act_ste | K_act_ste, diag_matrix(rep_vector(1.0, K_act_ste)));
//   target += inv_wishart_lpdf(cov_i_mono | K_mono, diag_matrix(rep_vector(1.0, K_mono)));
//   target += inv_wishart_lpdf(cov_i_act_ste | K_act_ste, diag_matrix(rep_vector(1.0, K_act_ste)));

  sigma_mono ~ cauchy(0, 5); // prior for sigma
  Lcorr_mono ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
  sigma_act_ste ~ cauchy(0, 5); // prior for sigma
  Lcorr_act_ste ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix

  target += normal_lpdf(gamma | 0, 1);
  target += inv_wishart_lpdf(cov_e | K, diag_matrix(rep_vector(1.0, K)));


  //likelihood
//   for (p in 1:P_mono) {
// 	target += multi_normal_lpdf(mu_mono[p] | eta_mono, cov_e_mono); 
//   }
  for (n in 1:N_mono) {
	target += multi_normal_lpdf(y_mono[n] | mu_mono[x_mono[n]], cov_i_mono);
  }
  
//   for (p in 1:P_act_ste) {
// 	target += multi_normal_lpdf(mu_act_ste[p] | eta_act_ste, cov_e_act_ste); 
//   }
  for (n in 1:N_act_ste) {
	target += multi_normal_lpdf(y_act_ste[n] | mu_act_ste[x_act_ste[n]], cov_i_act_ste);
  }

  //likelihood
  for (p in P_list) {
     target += multi_normal_lpdf(theta[p] | gamma, cov_e);
  }
  
}

generated quantities {
  matrix[N_mono, K_mono] Y_mono;     // monocyte signature loading matrix
  matrix[N_act_ste, K_act_ste] Y_act_ste;     // activated stellate signature loading matrix

  for (n in 1:N_mono) {
	Y_mono[n] = to_row_vector(multi_normal_rng(mu_mono[x_mono[n]], cov_i_mono));
  }
  for (n in 1:N_act_ste) {
	Y_act_ste[n] = to_row_vector(multi_normal_rng(mu_act_ste[x_act_ste[n]], cov_i_act_ste));
  }

} // The posterior predictive distribution