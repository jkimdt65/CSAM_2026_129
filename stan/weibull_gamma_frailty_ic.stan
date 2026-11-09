// -----------------------------------------------------------------
// weibull_gamma_frailty_ic.stan
//
// Replication of Goethals et al. (2009) / Gamma_Frailty_Interval_Censoring()
//   H(t | x, w) = w * lambda * t^gamma * exp(x*beta)
//   w_i ~ Gamma(1/theta, 1/theta)   (mean 1, variance theta)
//
// Key difference from the original R implementation:
//   The R code integrates out the gamma frailty ANALYTICALLY, which requires
//   summing 2^(#events in cluster) inclusion-exclusion terms (the `signs`
//   kronecker construction) and is hard-coded up to 4 events per cluster.
//   Stan instead samples the frailty w_i directly as a continuous latent
//   variable -- linear cost in cluster size, no combinatorial expansion,
//   and no restriction to the gamma family.
//
// Time scale: t is in QUARTER-YEARS (original days / 91.31), matching the
// R implementation, so estimates are directly comparable.
//
// censoring coding (cens_type):
//   2 = right censored   -> quarter infection-free through time L
//   4 = interval censored -> infection occurred in (L, U]
// -----------------------------------------------------------------

functions {
  // Weibull cumulative baseline hazard, safe at t = 0.
  //
  // Why this matters here: 68 quarters were already positive at the FIRST
  // inspection, so L = 0 for them (left-censored). Mathematically
  // H0(0) = lambda * 0^gamma = 0, and log_diff_exp(0, -HU) then reduces
  // exactly to log{1 - S(U)}, the correct left-censored contribution --
  // so no separate censoring branch is needed. But pow(0, gamma) has an
  // ill-defined derivative w.r.t. gamma in reverse-mode autodiff
  // (0 * log(0) -> NaN), which would poison the gradient. Returning a
  // hard 0 short-circuits the autodiff path.
  real weib_cumhaz(real t, real lambda, real gamma_shape) {
    if (t <= 0) return 0;
    return lambda * pow(t, gamma_shape);
  }
}

data {
  int<lower=1> N;                        // number of udder quarters
  int<lower=1> n_clusters;               // number of cows
  array[N] int<lower=1> cluster;         // cow index for each quarter
  array[N] int<lower=2, upper=4> cens_type;
  vector<lower=0>[N] L;                  // for type 4: last negative time
                                         // for type 2: censoring time
  vector<lower=0>[N] U;                  // for type 4: first positive time
  vector[N] x;                           // covariate: 1 if lactation > 1
}

parameters {
  real<lower=0> theta;                   // gamma frailty variance
  real<lower=0> lambda;                  // Weibull scale
  real<lower=0> gamma_shape;             // Weibull shape
  real beta;                             // log hazard ratio for x
  vector<lower=0>[n_clusters] w;         // frailties, sampled directly
}

model {
  // ---- priors ----
  // weakly informative; the data are highly informative here (N = 4784)
  theta       ~ normal(0, 5);
  lambda      ~ normal(0, 5);
  gamma_shape ~ normal(0, 5);
  beta        ~ normal(0, 5);

  // gamma frailty with mean 1, variance theta
  w ~ gamma(1 / theta, 1 / theta);

  // ---- likelihood ----
  for (i in 1:N) {
    real ex = exp(x[i] * beta);
    real HL = w[cluster[i]] * weib_cumhaz(L[i], lambda, gamma_shape) * ex;

    if (cens_type[i] == 2) {
      // right censored: log S(L) = -H(L)
      target += -HL;
    } else {
      // interval censored: log{ S(L) - S(U) }
      // when L = 0 this is log{1 - S(U)}, i.e. left-censored -- handled
      // by the same expression, no separate branch required
      real HU = w[cluster[i]] * weib_cumhaz(U[i], lambda, gamma_shape) * ex;
      target += log_diff_exp(-HL, -HU);
    }
  }
}

generated quantities {
  // Kendall's tau for the gamma frailty model -- the standard summary of
  // within-cow dependence, directly comparable to the published analysis
  real kendall_tau = theta / (theta + 2);
  real hazard_ratio = exp(beta);

  vector[N] log_lik;
  for (i in 1:N) {
    real ex = exp(x[i] * beta);
    real HL = w[cluster[i]] * weib_cumhaz(L[i], lambda, gamma_shape) * ex;
    if (cens_type[i] == 2) {
      log_lik[i] = -HL;
    } else {
      real HU = w[cluster[i]] * weib_cumhaz(U[i], lambda, gamma_shape) * ex;
      log_lik[i] = log_diff_exp(-HL, -HU);
    }
  }
}
