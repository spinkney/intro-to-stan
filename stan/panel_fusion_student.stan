data {
  int<lower=0> N;
  int<lower=0> N_unknown;
  vector<lower=0>[N] y_small;
  vector<lower=0>[N - N_unknown] y;
  vector[N - N_unknown] hours_big;
  vector[N - N_unknown] n_big;
  
  vector[N] n_small;
  vector[N] hours_small;
  
  array[N_unknown] int unknown;
  array[N - N_unknown] int known;
  
}
transformed data {
  vector[N - N_unknown] std_big;
  real mean_big = mean(hours_big);
  real mean_small = mean(hours_small);
  real sd_big = sd(hours_big);
  real sd_small = sd(hours_small);
  vector[2] sd_hours = [sd_big, sd_small]';
  vector[2] mean_hours = [mean_big, mean_small]';
  
  matrix[N, 2] std_small;
  
  real max_hrs;
  real min_hrs;
  
  std_big = ( hours_big - mean_big ) / sd_big;
  std_small[, 1] = ( n_small - mean(n_small) ) / sd(n_small);
  std_small[, 2] = ( hours_small - mean_small ) / sd_small;
  
  max_hrs = max(std_big);
  min_hrs = min(std_big);
}
parameters {
  vector<lower=0>[N_unknown] y_pred;
  vector<lower=min_hrs, upper=max_hrs>[N_unknown] x_pred;
  
  vector[3] beta_x;
  vector[2] alpha_x;
  real<lower=2> nu;
  
  vector[2] alpha;
  vector[3] beta;

  cholesky_factor_corr[2] L;
  vector<lower=0>[2] sigma;
  
  real<lower=0, upper=1> rho_xy;
}
transformed parameters {
  array[2] vector[N] x_all_raw;
  x_all_raw[1, unknown] = x_pred;
  x_all_raw[1, known] = std_big;
  x_all_raw[2] = std_small[, 2];
}
model {
  // priors
  beta_x ~ std_normal();
  alpha_x ~ normal(0, 4);
 
  alpha ~ normal(0, 4);
  beta ~ std_normal();
   
  sigma ~ exponential(0.25);

  L ~ lkj_corr_cholesky(1);
  nu ~ exponential(1);
  rho_xy ~ beta(2, 2);
  
  // likelihood
  array[N] vector[2] y_all;
  array[N] vector[2] x_all;
  vector[N] y_big;
  y_big[known] = y;
  y_big[unknown] = y_pred;
   
  array[N] vector[2] mu;
  array[N] vector[2] mu_x;
  for (i in 1:N) {
    x_all[i, 1] = x_all_raw[1, i];
    x_all[i, 2] = x_all_raw[2, i];
    y_all[i, 1] = y_big[i];
    y_all[i, 2] = y_small[i];
     
    mu[i, 1] = alpha[1] + x_all[i, 1] * beta[1] + x_all[i, 2] * beta[2];
    mu[i, 2] = alpha[2] + x_all[i, 2] * beta[3];
    
    mu_x[i, 1] = alpha_x[1] + y_all[i, 1] * beta_x[1] + y_all[i, 2] * beta_x[2];
    mu_x[i, 2] = alpha_x[2] + y_all[i, 2] * beta_x[3];
  }
  x_all ~ multi_normal_cholesky(mu_x, rho_xy * L);
  y_all ~ multi_student_t_cholesky(nu, mu, rho_xy * diag_pre_multiply(sigma, L));
}
generated quantities {
   vector[N_unknown] x_out = x_pred * sd_big + mean_big;
   vector[N_unknown] y_out = y_pred * 400;
   
   matrix[N, 2] x_rep;
   matrix[N, 2] y_rep;
   
   {
    array[N] vector[2] y_all;
    vector[N] y_big;
    y_big[known] = y;
    y_big[unknown] = y_pred;
   
    array[N] vector[2] mu;
     array[N] vector[2] mu_x;
    for (i in 1:N) {
      y_all[i, 1] = y_big[i];
      y_all[i, 2] = y_small[i];
      mu[i, 1] = alpha[1] + x_all_raw[1, i] * beta[1] + x_all_raw[2, i] * beta[2];
      mu[i, 2] = alpha[2] + x_all_raw[2, i] * beta[3];
      mu_x[i, 1] = alpha_x[1] + y_all[i, 1] * beta_x[1] + y_all[i, 2] * beta_x[2];
      mu_x[i, 2] = alpha_x[2] + y_all[i, 2] * beta_x[3];
      
      x_rep[i] = (multi_normal_cholesky_rng(mu_x[i], rho_xy * L) .* sd_hours + mean_hours)';
       while (x_rep[i, 1] < 0 || x_rep[i, 2] < 0 )
          x_rep[i] = (multi_normal_cholesky_rng(mu_x[i], rho_xy * L) .* sd_hours + mean_hours)';

      y_rep[i] = multi_student_t_cholesky_rng(nu, mu[i], rho_xy * diag_pre_multiply(sigma, L))' .* [400, 10];
      while (y_rep[i, 1] < 0 || y_rep[i, 2] < 0) 
         y_rep[i] = multi_student_t_cholesky_rng(nu, mu[i], rho_xy * diag_pre_multiply(sigma, L))' .* [400, 10];
    }
  }
}
