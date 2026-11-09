// -----------------------------------------------------------------
// weibull_gamma_frailty_X.stan
// Weibull baseline + cow-level gamma frailty with a covariate matrix X
// (Table 4, Weibull column; within- vs between-cluster comparison).
//
// Same model as weibull_gamma_frailty_ic.stan (Table 2) with the single
// covariate x replaced by a matrix X, so that the three baselines of
// Table 4 receive the same covariates.
//
// Identifiability: lambda carries the level of the baseline and E[w] = 1
// fixes the level of the frailties, so no alpha/xi split is needed here.
// X must not contain an intercept column (confounded with lambda).
// -----------------------------------------------------------------

functions {
  real weib_cumhaz(real t, real lambda, real gamma_shape) {
    if (t <= 0) return 0;
    return lambda * pow(t, gamma_shape);
  }
}

data {
  int<lower=1> N;
  int<lower=0> P;
  matrix[N, P] X;                        // no intercept column

  int<lower=1> n_clusters;
  array[N] int<lower=1> cluster;
  array[N] int<lower=2, upper=4> cens_type;
  vector<lower=0>[N] L;
  vector<lower=0>[N] U;
}

parameters {
  real<lower=0> theta;                   // gamma frailty variance
  real<lower=0> lambda;
  real<lower=0> gamma_shape;
  vector[P] beta;
  vector<lower=0>[n_clusters] w;
}

model {
  theta       ~ normal(0, 5);
  lambda      ~ normal(0, 5);
  gamma_shape ~ normal(0, 5);
  beta        ~ normal(0, 5);

  w ~ gamma(1 / theta, 1 / theta);

  for (i in 1:N) {
    real ex = exp(X[i] * beta);
    real HL = w[cluster[i]] * weib_cumhaz(L[i], lambda, gamma_shape) * ex;

    if (cens_type[i] == 2) {
      target += -HL;
    } else {
      real HU = w[cluster[i]] * weib_cumhaz(U[i], lambda, gamma_shape) * ex;
      target += log_diff_exp(-HL, -HU);
    }
  }
}

generated quantities {
  vector[P] hazard_ratio = exp(beta);
  real kendall_tau = theta / (theta + 2);
  real var_log_frailty = trigamma(1 / theta);
  real mean_log_w = mean(log(w));

  vector[N] log_lik;
  for (i in 1:N) {
    real ex = exp(X[i] * beta);
    real HL = w[cluster[i]] * weib_cumhaz(L[i], lambda, gamma_shape) * ex;
    if (cens_type[i] == 2) {
      log_lik[i] = -HL;
    } else {
      real HU = w[cluster[i]] * weib_cumhaz(U[i], lambda, gamma_shape) * ex;
      log_lik[i] = log_diff_exp(-HL, -HU);
    }
  }
}
