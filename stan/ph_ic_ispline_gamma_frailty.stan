// -----------------------------------------------------------------
// ph_ic_ispline_gamma_frailty.stan
// I-spline baseline + cow-level gamma frailty (Table 4, Supplementary
// Table S1).
//   H0(t)      = sum_k gamma_k * I_k(t),   gamma_k > 0
//   H_ij(t|w)  = w_i * H0(t) * exp(x_ij' beta)
//   w_i ~ Gamma(1/theta, 1/theta)
//
// Monotonicity holds by construction (gamma_k > 0, I_k non-decreasing),
// and since I_k(0) = 0, rows with L = 0 give H0(L) = 0 and the
// interval-censored contribution reduces to the left-censored one.
//
// The M-spline basis for h0(t) is needed only for exact event times
// (cens_type == 1), which the udder quarter data do not contain, so it
// is not passed.
//
// Identifiability, as in ph_ic_pwe_gamma_frailty.stan:
//   log gamma_k = alpha + xi[k],   sum(xi) = 0,
// with alpha carrying the level and the random walk acting on xi only.
//
// The basis matrices are computed in R with splines2 and passed as data.
// -----------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=0> P;
  matrix[N, P] X;

  int<lower=1> n_clusters;
  array[N] int<lower=1> cluster;

  array[N] int<lower=1, upper=4> cens_type;
  vector<lower=0>[N] entry;              // all 0 for the udder quarter data

  int<lower=2> nbasis;
  matrix[N, nbasis] IL;                  // I-spline basis at L
  matrix[N, nbasis] IU;                  // I-spline basis at U
  matrix[N, nbasis] Ient;                // at entry (I(0) = 0 when entry = 0)
}

parameters {
  vector[P] beta;

  real alpha;                            // level of the log coefficients
  vector[nbasis - 1] xi_free;            // shape
  real<lower=0> smooth_sd;

  real<lower=0> theta;                   // gamma frailty variance
  vector<lower=0>[n_clusters] w;
}

transformed parameters {
  vector[nbasis] xi;
  vector[nbasis] gamma_coef;

  xi[1:(nbasis - 1)] = xi_free;
  xi[nbasis] = -sum(xi_free);            // sum(xi) = 0
  gamma_coef = exp(alpha + xi);          // > 0
}

model {
  beta  ~ normal(0, 5);
  alpha ~ normal(0, 5);

  smooth_sd ~ normal(0, 1);
  for (k in 2:nbasis) {
    xi[k] ~ normal(xi[k - 1], smooth_sd);   // random walk on the shape only
  }

  theta ~ normal(0, 5);
  w ~ gamma(1 / theta, 1 / theta);

  for (i in 1:N) {
    real eta  = X[i] * beta;
    real mult = w[cluster[i]] * exp(eta);
    real logS_entry = 0;

    if (entry[i] > 0) {
      logS_entry = -mult * dot_product(Ient[i], gamma_coef);
    }

    if (cens_type[i] == 2) {
      target += -mult * dot_product(IL[i], gamma_coef) - logS_entry;

    } else if (cens_type[i] == 3) {
      target += log1m_exp(-mult * dot_product(IU[i], gamma_coef)) - logS_entry;

    } else {
      // interval-censored; L = 0 gives a zero IL row, H_L = 0, i.e. left-censored
      real H_L = dot_product(IL[i], gamma_coef);
      real H_U = dot_product(IU[i], gamma_coef);
      target += log_diff_exp(-mult * H_L, -mult * H_U) - logS_entry;
    }
  }
}

generated quantities {
  vector[P] hazard_ratio = exp(beta);
  real kendall_tau = theta / (theta + 2);
  real var_log_frailty = trigamma(1 / theta);
  real mean_log_w = mean(log(w));        // for the correlation with alpha

  vector[N] log_lik;
  for (i in 1:N) {
    real eta  = X[i] * beta;
    real mult = w[cluster[i]] * exp(eta);
    real logS_entry = 0;

    if (entry[i] > 0) {
      logS_entry = -mult * dot_product(Ient[i], gamma_coef);
    }

    if (cens_type[i] == 2) {
      log_lik[i] = -mult * dot_product(IL[i], gamma_coef) - logS_entry;

    } else if (cens_type[i] == 3) {
      log_lik[i] = log1m_exp(-mult * dot_product(IU[i], gamma_coef)) - logS_entry;

    } else {
      real H_L = dot_product(IL[i], gamma_coef);
      real H_U = dot_product(IU[i], gamma_coef);
      log_lik[i] = log_diff_exp(-mult * H_L, -mult * H_U) - logS_entry;
    }
  }
}
