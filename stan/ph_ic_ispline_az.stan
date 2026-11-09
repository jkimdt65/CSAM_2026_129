// ph_ic_ispline_az.stan
// I-spline baseline, PH model, interval-censored data, no frailty.
// Log coefficients written as alpha + xi with sum(xi) = 0 (Section 2.2);
// RW1 smoothing on xi alone. Basis built in R.
// Prior scales passed as data.
//
// cens_type: 1 = exact, 2 = right, 3 = left, 4 = interval

data {
  int<lower=1> N;
  int<lower=0> P;
  matrix[N, P] X;
  array[N] int<lower=1, upper=4> cens_type;
  vector<lower=0>[N] entry;
  int<lower=2> nbasis;
  matrix[N, nbasis] IL;
  matrix[N, nbasis] IU;
  matrix[N, nbasis] Ient;
  matrix[N, nbasis] ML;
  real<lower=0> beta_sd;
  real<lower=0> smooth_scale;
}

parameters {
  vector[P] beta;
  real alpha;
  vector[nbasis - 1] xi_free;
  real<lower=0> smooth_sd;
}

transformed parameters {
  vector[nbasis] xi = append_row(xi_free, -sum(xi_free));
  vector[nbasis] gamma = exp(alpha + xi);               // positive by construction
}

model {
  beta ~ normal(0, beta_sd);
  alpha ~ normal(0, 5);
  smooth_sd ~ normal(0, smooth_scale);
  for (k in 2:nbasis) xi[k] ~ normal(xi[k - 1], smooth_sd);

  for (i in 1:N) {
    real eta = X[i] * beta;
    real ex = exp(eta);
    real logS_entry = -ex * dot_product(Ient[i], gamma);
    if (cens_type[i] == 1) {
      target += log(dot_product(ML[i], gamma)) + eta
                - ex * dot_product(IL[i], gamma) - logS_entry;
    } else if (cens_type[i] == 2) {
      target += -ex * dot_product(IL[i], gamma) - logS_entry;
    } else if (cens_type[i] == 3) {
      target += log1m_exp(-ex * dot_product(IU[i], gamma)) - logS_entry;
    } else {
      target += log_diff_exp(-ex * dot_product(IL[i], gamma),
                             -ex * dot_product(IU[i], gamma)) - logS_entry;
    }
  }
}
