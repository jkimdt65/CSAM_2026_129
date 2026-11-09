// ph_ic_pwe_az.stan
// Piecewise exponential baseline, PH model, interval-censored data, no frailty.
// The log hazards are written as alpha + xi with sum(xi) = 0 (Section 2.2). alpha carries the level,
// xi the shape, and the RW1 smoothing prior applies to xi alone.
// Prior scales are passed as data so that the same program serves the main
// fit and the prior sensitivity analysis.
//
// cens_type: 1 = exact, 2 = right, 3 = left, 4 = interval

functions {
  real pw_cumhaz(real t, vector cuts, vector log_lambda) {
    int K = num_elements(log_lambda);
    real H = 0;
    if (t <= 0) return 0;
    for (k in 1:K) {
      if (t > cuts[k]) H += exp(log_lambda[k]) * (fmin(t, cuts[k + 1]) - cuts[k]);
    }
    return H;
  }
  real pw_loghaz(real t, vector cuts, vector log_lambda) {
    int K = num_elements(log_lambda);
    for (k in 1:K) if (t > cuts[k] && t <= cuts[k + 1]) return log_lambda[k];
    return log_lambda[K];
  }
}

data {
  int<lower=1> N;
  int<lower=0> P;
  matrix[N, P] X;
  array[N] int<lower=1, upper=4> cens_type;
  vector<lower=0>[N] L;
  vector<lower=0>[N] U;
  vector<lower=0>[N] entry;
  int<lower=2> K;
  vector[K + 1] cuts;
  real<lower=0> beta_sd;        // prior sd of beta (main analysis: 5)
  real<lower=0> smooth_scale;   // smooth_sd ~ half-N(0, smooth_scale) (main: 1)
}

parameters {
  vector[P] beta;
  real alpha;
  vector[K - 1] xi_free;
  real<lower=0> smooth_sd;
}

transformed parameters {
  vector[K] xi = append_row(xi_free, -sum(xi_free));   // sums to zero
  vector[K] log_lambda = alpha + xi;                   // level + shape
}

model {
  beta ~ normal(0, beta_sd);
  alpha ~ normal(0, 5);
  smooth_sd ~ normal(0, smooth_scale);
  for (k in 2:K) xi[k] ~ normal(xi[k - 1], smooth_sd);  // RW1 on shape only

  for (i in 1:N) {
    real eta = X[i] * beta;
    real ex = exp(eta);
    real logS_entry = entry[i] > 0 ? -ex * pw_cumhaz(entry[i], cuts, log_lambda) : 0;
    if (cens_type[i] == 1) {
      target += pw_loghaz(L[i], cuts, log_lambda) + eta
                - ex * pw_cumhaz(L[i], cuts, log_lambda) - logS_entry;
    } else if (cens_type[i] == 2) {
      target += -ex * pw_cumhaz(L[i], cuts, log_lambda) - logS_entry;
    } else if (cens_type[i] == 3) {
      target += log1m_exp(-ex * pw_cumhaz(U[i], cuts, log_lambda)) - logS_entry;
    } else {
      target += log_diff_exp(-ex * pw_cumhaz(L[i], cuts, log_lambda),
                             -ex * pw_cumhaz(U[i], cuts, log_lambda)) - logS_entry;
    }
  }
}
