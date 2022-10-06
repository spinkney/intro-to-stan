functions {
  real inv_Phi_log_lambda(real log_p) {
    vector[8] a
        = [3.3871328727963666080e+00, 1.3314166789178437745e+02,
           1.9715909503065514427e+03, 1.3731693765509461125e+04,
           4.5921953931549871457e+04, 6.7265770927008700853e+04,
           3.3430575583588128105e+04, 2.5090809287301226727e+03]';
    vector[7] b
        = [4.2313330701600911252e+01, 6.8718700749205790830e+02,
           5.3941960214247511077e+03, 2.1213794301586595867e+04,
           3.9307895800092710610e+04, 2.8729085735721942674e+04,
           5.2264952788528545610e+03]';
    vector[8] c
        = [1.42343711074968357734e+00, 4.63033784615654529590e+00,
           5.76949722146069140550e+00, 3.64784832476320460504e+00,
           1.27045825245236838258e+00, 2.41780725177450611770e-01,
           2.27238449892691845833e-02, 7.74545014278341407640e-04]';
    vector[7] d
        = [2.05319162663775882187e+00, 1.67638483018380384940e+00,
           6.89767334985100004550e-01, 1.48103976427480074590e-01,
           1.51986665636164571966e-02, 5.47593808499534494600e-04,
           1.05075007164441684324e-09]';
    vector[8] e
        = [6.65790464350110377720e+00, 5.46378491116411436990e+00,
           1.78482653991729133580e+00, 2.96560571828504891230e-01,
           2.65321895265761230930e-02, 1.24266094738807843860e-03,
           2.71155556874348757815e-05, 2.01033439929228813265e-07]';
    vector[7] f
        = [5.99832206555887937690e-01, 1.36929880922735805310e-01,
           1.48753612908506148525e-02, 7.86869131145613259100e-04,
           1.84631831751005468180e-05, 1.42151175831644588870e-07,
           2.04426310338993978564e-15]';
    real log_q = log_p <= log(0.5) ? log_diff_exp(log(1), log_sum_exp(log_p, log(0.5))) : log_diff_exp(log_p, log(0.5));
    int log_q_sign = log_p <= log(0.5) ? -1 : 1;
    real val;

    if (log_p == log(1)) {
      return positive_infinity();
    }

    if (log_q <= log(.425)) {
      real log_r;
      real log_rtn;
      real log_agg_a;
      real log_agg_b;
      vector[8] log_a = log(a);
      vector[7] log_b = log(b);
      log_r = log_diff_exp(log(.180625), 2 * log_q);
      log_agg_a = log_sum_exp(log_a[8] + log_r, log_a[7]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[6]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[5]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[4]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[3]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[2]);
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[1]);

      log_agg_b = log_sum_exp(log_b[7] + log_r, log_b[6]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[5]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[4]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[3]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[2]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[1]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, 0);

      log_rtn = log_q + log_agg_a - log_agg_b;
      return log_q_sign * exp(log_rtn);
    } else {
      real log_r = log_q_sign == -1 ? log_p : log_diff_exp(log(1), log_p);

      if (is_inf(log_r)) {
        return 0;
      }

      log_r = log(sqrt(-log_r));

      if (log_r <= log(5.0)) {
        vector[8] log_c = log(c);
        vector[7] log_d = log(d);
        real log_agg_c;
        real log_agg_d;
        log_r = log_diff_exp(log_r, log(1.6));

        log_agg_c = log_sum_exp(log_c[8] + log_r, log_c[7]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[6]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[5]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[4]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[3]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[2]);
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[1]);

        log_agg_d = log_sum_exp(log_d[7] + log_r, log_d[6]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[5]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[4]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[3]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[2]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[1]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, 0);

        val = exp(log_agg_c - log_agg_d);
      } else {
        vector[8] log_e = log(e);
        vector[7] log_f = log(f);
        real log_agg_e;
        real log_agg_f;
        log_r = log_diff_exp(log_r, log(5));

        log_agg_e = log_sum_exp(log_e[8] + log_r, log_e[7]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[6]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[5]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[4]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[3]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[2]);
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[1]);

        log_agg_f = log_sum_exp(log_f[7] + log_r, log_f[6]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[5]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[4]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[3]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[2]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[1]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, 0);

        val = exp(log_agg_e - log_agg_f);
      }
      if (log_q_sign == -1)
        return -val;
    }
    return val;
  }

  real inv_Phi_log_fun(real log_p) {
    real log_BIGINT = log(2000000000);
    return log_p >= log(0.9999) ? -inv_Phi_log_lambda(
               log_diff_exp(log_BIGINT, log_BIGINT + log_p) - log_BIGINT)
                       : inv_Phi_log_lambda(log_p);
  }
   real binormal_lccdf(real z, real rho) {
    return log(0.75 - asin(rho) / (2 * pi()));
  }
  
real multi_normal_cholesky_copula_lpdf(matrix U, matrix L) {
  int N = rows(U);
  int J = cols(U);
  matrix[J, J] Gammainv = chol2inv(L);
  return -N * sum(log(diagonal(L))) - 0.5 * sum(add_diag(Gammainv, -1.0) .* crossprod(U));
}
}
data {
  int<lower=0> N;
  int<lower=0> N_unknown;
  vector<lower=0>[N] y_noisy;
  vector<lower=0>[N - N_unknown] y;
  vector[N - N_unknown] hours_inscape;
  vector[N - N_unknown] n_inscape;
  
  vector[N] n_thp;
  vector[N] hours_thp;
  
  array[N_unknown] int unknown;
  array[N - N_unknown] int known;
  
}
transformed data {
  vector[N - N_unknown] std_inscape;
  real mean_inscape = mean(hours_inscape);
  real sd_inscape = sd(hours_inscape);
  
  matrix[N, 2] std_thp;
  
  real max_hrs;
  real min_hrs;
  
  std_inscape = ( hours_inscape - mean_inscape ) / sd_inscape;
  std_thp[, 1] = ( n_thp - mean(n_thp) ) / sd(n_thp);
  std_thp[, 2] = ( hours_thp - mean(hours_thp) ) / sd(hours_thp);
  
  max_hrs = max(std_inscape);
  min_hrs = min(std_inscape);
}
parameters {
  vector<lower=0>[N_unknown] y_pred;
  vector<lower=min_hrs, upper=max_hrs>[N_unknown] x_pred;
  
  vector[2] beta_x;
  real<lower=0> sigma_x;
  real alpha_x;
  
  vector[2] alpha;
  vector[2] beta;
  
  cholesky_factor_corr[2] L;
  vector<lower=0>[2] sigma;
}
transformed parameters {
  vector[N] x_all;
  x_all[unknown] = x_pred;
  x_all[known] = std_inscape;
}
model {
  // priors
   beta_x ~ std_normal();
   alpha_x ~ std_normal();
   sigma_x ~ exponential(1);
   
   alpha ~ normal(0, 4);
   beta ~ std_normal();
   
   sigma ~ exponential(1);
   L ~ lkj_corr_cholesky(1);
  
  // likelihood
   x_all ~ normal_id_glm(std_thp, alpha_x, beta_x, sigma_x);
   
   vector[N] y_big;
   y_big[known] = y;
   y_big[unknown] = y_pred;
   
   matrix[N, 2] y_normal; 
   array[2] vector[N] mu;
   mu[1] = alpha[1] + x_all * beta[1];
   mu[2] = alpha[2] + std_thp[, 2] * beta[2];
   
   y_big ~ normal(mu[1] , sigma[1]);
   y_noisy ~ normal(mu[2] , sigma[2]);
   for (i in 1:N) {
     y_normal[i, 1] = inv_Phi_log_fun(normal_lcdf(y_big[i] | mu[1, i], sigma[1]));
     y_normal[i, 2] = inv_Phi_log_fun(normal_lcdf(y_noisy[i] | mu[2, i], sigma[2]));
   }
 
  y_normal ~ multi_normal_cholesky_copula(L);
  target += -N * binormal_lccdf(0.0 | L[2, 1]); 

}
generated quantities {
   vector[N_unknown] x_out = x_pred * sd_inscape + mean_inscape;
   vector[N_unknown] y_out = y_pred * 400;
}
