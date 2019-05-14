// generated with brms 2.8.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  int<lower=1> N_Y1;  // number of observations
  vector[N_Y1] Y_Y1;  // response variable
  int<lower=1> K_Y1;  // number of population-level effects
  matrix[N_Y1, K_Y1] X_Y1;  // population-level design matrix
  int<lower=1> N_Y2;  // number of observations
  vector[N_Y2] Y_Y2;  // response variable
  int<lower=1> K_Y2;  // number of population-level effects
  matrix[N_Y2, K_Y2] X_Y2;  // population-level design matrix
  int<lower=1> nresp;  // number of responses
  int nrescor;  // number of residual correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_Y1 = K_Y1 - 1;
  matrix[N_Y1, K_Y1 - 1] Xc_Y1;  // centered version of X_Y1
  vector[K_Y1 - 1] means_X_Y1;  // column means of X_Y1 before centering
  int Kc_Y2 = K_Y2 - 1;
  matrix[N_Y2, K_Y2 - 1] Xc_Y2;  // centered version of X_Y2
  vector[K_Y2 - 1] means_X_Y2;  // column means of X_Y2 before centering
  vector[nresp] Y[N];  // response matrix
  for (i in 2:K_Y1) {
    means_X_Y1[i - 1] = mean(X_Y1[, i]);
    Xc_Y1[, i - 1] = X_Y1[, i] - means_X_Y1[i - 1];
  }
  for (i in 2:K_Y2) {
    means_X_Y2[i - 1] = mean(X_Y2[, i]);
    Xc_Y2[, i - 1] = X_Y2[, i] - means_X_Y2[i - 1];
  }
  for (n in 1:N) {
    Y[n] = [Y_Y1[n], Y_Y2[n]]';
  }
}
parameters {
  vector[Kc_Y1] b_Y1;  // population-level effects
  real temp_Y1_Intercept;  // temporary intercept
  real<lower=0> sigma_Y1;  // residual SD
  vector[Kc_Y2] b_Y2;  // population-level effects
  real temp_Y2_Intercept;  // temporary intercept
  real<lower=0> sigma_Y2;  // residual SD
  // parameters for multivariate linear models
  cholesky_factor_corr[nresp] Lrescor;
}
transformed parameters {
}
model {
  vector[N_Y1] mu_Y1 = temp_Y1_Intercept + Xc_Y1 * b_Y1;
  vector[N_Y2] mu_Y2 = temp_Y2_Intercept + Xc_Y2 * b_Y2;
  // multivariate linear predictor matrix
  vector[nresp] Mu[N];
  vector[nresp] sigma = [sigma_Y1, sigma_Y2]';
  // cholesky factor of residual covariance matrix
  matrix[nresp, nresp] LSigma = diag_pre_multiply(sigma, Lrescor);
  for (n in 1:N) {
    Mu[n] = [mu_Y1[n], mu_Y2[n]]';
  }
  // priors including all constants
  target += student_t_lpdf(temp_Y1_Intercept | 3, 1, 10);
  target += student_t_lpdf(sigma_Y1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(temp_Y2_Intercept | 3, 5, 10);
  target += student_t_lpdf(sigma_Y2 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += lkj_corr_cholesky_lpdf(Lrescor | 1);
  // likelihood including all constants
  if (!prior_only) {
    target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Y1_Intercept = temp_Y1_Intercept - dot_product(means_X_Y1, b_Y1);
  // actual population-level intercept
  real b_Y2_Intercept = temp_Y2_Intercept - dot_product(means_X_Y2, b_Y2);
  corr_matrix[nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor);
  vector<lower=-1,upper=1>[nrescor] rescor;
  // extract upper diagonal of correlation matrix
  for (k in 1:nresp) {
    for (j in 1:(k - 1)) {
      rescor[choose(k - 1, 2) + j] = Rescor[j, k];
    }
  }
}
