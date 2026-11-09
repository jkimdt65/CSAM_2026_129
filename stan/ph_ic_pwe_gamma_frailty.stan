// -----------------------------------------------------------------
// ph_ic_pwe_gamma_frailty.stan
// Piecewise exponential baseline + cow-level gamma frailty (Table 4,
// Supplementary Table S1).
//   h_ij(t | w) = h0(t) * w_i * exp(x_ij' beta)
//   w_i ~ Gamma(1/theta, 1/theta)      (mean 1, variance theta)
//
// Here theta denotes the frailty variance; the local variable holding
// w_i * exp(eta) is called mult.
//
// Identifiability (Section 2.2). Dividing h0 by c and multiplying every
// w_i by c leaves the likelihood unchanged. A random-walk prior acts on
// differences and hardly constrains the level, so the log hazards are
// written as
//   log_lambda[k] = alpha + xi[k],   sum(xi) = 0,
// with alpha carrying the level and the random walk acting on xi only.
// sum(xi) = 0 is imposed by xi[K] = -sum(xi[1:(K-1)]), a linear map with
// constant Jacobian.
// -----------------------------------------------------------------

functions {
  real pw_cumhaz(real t, vector cuts, vector log_lambda) {
    int K = num_elements(log_lambda);
    real H = 0;
    if (t <= 0) return 0;
    for (k in 1:K) {
      real lo = cuts[k];
      real hi = cuts[k + 1];
      if (t > lo) {
        real seg = fmin(t, hi) - lo;
        H += exp(log_lambda[k]) * seg;
      }
    }
    return H;
  }

  real pw_loghaz(real t, vector cuts, vector log_lambda) {
    int K = num_elements(log_lambda);
    for (k in 1:K) {
      if (t > cuts[k] && t <= cuts[k + 1]) {
        return log_lambda[k];
      }
    }
    return log_lambda[K];
  }
}

data {
  int<lower=1> N;
  int<lower=0> P;
  matrix[N, P] X;                     // no intercept column (alpha plays that role)

  int<lower=1> n_clusters;
  array[N] int<lower=1> cluster;

  array[N] int<lower=1, upper=4> cens_type;
  vector<lower=0>[N] L;
  vector<lower=0>[N] U;
  vector<lower=0>[N] entry;           // all 0 for the udder quarter data

  int<lower=2> K;
  vector[K + 1] cuts;
}

parameters {
  vector[P] beta;

  real alpha;                         // level of the log baseline hazard
  vector[K - 1] xi_free;              // shape
  real<lower=0> smooth_sd;

  real<lower=0> theta;                // gamma frailty variance
  vector<lower=0>[n_clusters] w;      // frailties, sampled directly
}

transformed parameters {
  vector[K] xi;
  vector[K] log_lambda;

  xi[1:(K - 1)] = xi_free;
  xi[K] = -sum(xi_free);
  log_lambda = alpha + xi;
}

model {
  // ---- priors ----
  beta  ~ normal(0, 5);

  // alpha: level of the log baseline hazard. With time in quarter-years
  // the Weibull fit gives h0(1) of about 0.24 (log about -1.4), so
  // normal(0, 5) is wide.
  alpha ~ normal(0, 5);

  smooth_sd ~ normal(0, 1);
  for (k in 2:K) {
    xi[k] ~ normal(xi[k - 1], smooth_sd);   // random walk on the shape only
  }

  theta ~ normal(0, 5);
  w ~ gamma(1 / theta, 1 / theta);

  // ---- likelihood ----
  for (i in 1:N) {
    real eta  = X[i] * beta;
    real mult = w[cluster[i]] * exp(eta);   // the frailty enters here
    real logS_entry = 0;

    if (entry[i] > 0) {
      logS_entry = -mult * pw_cumhaz(entry[i], cuts, log_lambda);
    }

    if (cens_type[i] == 1) {
      real H_t   = pw_cumhaz(L[i], cuts, log_lambda);
      real logh0 = pw_loghaz(L[i], cuts, log_lambda);
      target += logh0 + eta + log(w[cluster[i]]) - mult * H_t - logS_entry;

    } else if (cens_type[i] == 2) {
      target += -mult * pw_cumhaz(L[i], cuts, log_lambda) - logS_entry;

    } else if (cens_type[i] == 3) {
      real logS_U = -mult * pw_cumhaz(U[i], cuts, log_lambda);
      target += log1m_exp(logS_U) - logS_entry;

    } else {
      real logS_L = -mult * pw_cumhaz(L[i], cuts, log_lambda);
      real logS_U = -mult * pw_cumhaz(U[i], cuts, log_lambda);
      target += log_diff_exp(logS_L, logS_U) - logS_entry;
    }
  }
}

generated quantities {
  vector[P] hazard_ratio = exp(beta);
  real kendall_tau = theta / (theta + 2);
  real var_log_frailty = trigamma(1 / theta);   // comparable to sigma_b^2 of the log-normal

  // A posterior correlation of alpha and mean_log_w near -1 would signal
  // that the level of the baseline and of the frailties is not separated.
  real mean_log_w = mean(log(w));

  vector[N] log_lik;
  for (i in 1:N) {
    real eta  = X[i] * beta;
    real mult = w[cluster[i]] * exp(eta);
    real logS_entry = 0;

    if (entry[i] > 0) {
      logS_entry = -mult * pw_cumhaz(entry[i], cuts, log_lambda);
    }

    if (cens_type[i] == 1) {
      real H_t   = pw_cumhaz(L[i], cuts, log_lambda);
      real logh0 = pw_loghaz(L[i], cuts, log_lambda);
      log_lik[i] = logh0 + eta + log(w[cluster[i]]) - mult * H_t - logS_entry;

    } else if (cens_type[i] == 2) {
      log_lik[i] = -mult * pw_cumhaz(L[i], cuts, log_lambda) - logS_entry;

    } else if (cens_type[i] == 3) {
      real logS_U = -mult * pw_cumhaz(U[i], cuts, log_lambda);
      log_lik[i] = log1m_exp(logS_U) - logS_entry;

    } else {
      real logS_L = -mult * pw_cumhaz(L[i], cuts, log_lambda);
      real logS_U = -mult * pw_cumhaz(U[i], cuts, log_lambda);
      log_lik[i] = log_diff_exp(logS_L, logS_U) - logS_entry;
    }
  }
}
