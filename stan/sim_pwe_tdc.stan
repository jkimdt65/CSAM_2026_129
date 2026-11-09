// sim_pwe_tdc.stan
// Simulation of Section 4: piecewise exponential baseline + time-varying
// covariate + delayed entry + (optional) log-normal frailty.
//
// Non-centred log-normal frailty:
//   log w_i = sigma_b * z_i,  z_i ~ N(0, 1)
// sigma_b and z are independent a priori, so no funnel arises.
// use_frailty = 0 fixes w = 1.
//
// Time-varying covariate: x2 switches from 0 to 1 at vchg (exogenous).
//   H(t) = w * exp(X*beta) * [ H0(min(t,v)) + (H0(t)-H0(v))_+ * exp(beta_tdc) ]
// Delayed entry: the cumulative hazard at entry is subtracted.

functions {
  // piecewise exponential cumulative hazard (Section 2.2)
  real pwe_cumhaz(real t, vector cuts, vector eta) {
    real H = 0;
    if (t <= 0) return 0;
    for (k in 1:num_elements(eta))
      if (t > cuts[k])
        H += exp(eta[k]) * (fmin(t, cuts[k + 1]) - cuts[k]);
    return H;
  }

  // subject-level cumulative hazard with the time-varying covariate
  real cumhaz_tdc(real t, real v, real lin, real beta_tdc,
                  vector cuts, vector eta) {
    real Hpre;
    real Hpost;
    if (t <= 0) return 0;
    Hpre  = pwe_cumhaz(fmin(t, v), cuts, eta);
    Hpost = t > v ? (pwe_cumhaz(t, cuts, eta) - pwe_cumhaz(v, cuts, eta)) : 0;
    return exp(lin) * (Hpre + Hpost * exp(beta_tdc));
  }
}

data {
  int<lower=1> N;
  int<lower=0> P;
  matrix[N, P] X;
  array[N] int<lower=2, upper=4> cens_type;   // 2 right, 3 left, 4 interval
  vector<lower=0>[N] L;
  vector[N] U;                                 // 1e4 for right censoring
  vector<lower=0>[N] entry;                    // left truncation time (0 if none)
  vector<lower=0>[N] vchg;                     // covariate change point

  int<lower=2> K;
  vector[K + 1] cuts;

  int<lower=1> n_clusters;
  array[N] int<lower=1> cluster;
  int<lower=0, upper=1> use_frailty;
}

parameters {
  vector[P] beta;
  real beta_tdc;
  real alpha;
  vector[K - 1] xi_free;
  real<lower=0> smooth_sd;

  array[use_frailty ? 1 : 0] real<lower=0> sigma_b;
  vector[use_frailty ? n_clusters : 0] z;
}

transformed parameters {
  vector[K] xi  = append_row(xi_free, -sum(xi_free));
  vector[K] eta = alpha + xi;

  vector[use_frailty ? n_clusters : 0] log_w;
  if (use_frailty) log_w = sigma_b[1] * z;
}

model {
  beta      ~ normal(0, 5);
  beta_tdc  ~ normal(0, 5);
  alpha     ~ normal(0, 5);
  smooth_sd ~ normal(0, 1);
  for (k in 2:K) xi[k] ~ normal(xi[k - 1], smooth_sd);

  if (use_frailty) {
    sigma_b[1] ~ normal(0, 1);
    z          ~ std_normal();
  }

  for (i in 1:N) {
    real wi  = use_frailty ? exp(log_w[cluster[i]]) : 1.0;
    real lin = P > 0 ? dot_product(X[i], beta) : 0.0;
    real HL  = wi * cumhaz_tdc(L[i], vchg[i], lin, beta_tdc, cuts, eta);
    real Hent = entry[i] > 0
              ? wi * cumhaz_tdc(entry[i], vchg[i], lin, beta_tdc, cuts, eta)
              : 0.0;

    if (cens_type[i] == 2) {
      target += -HL + Hent;                       // S(L)/S(A)
    } else {
      real HU = wi * cumhaz_tdc(U[i], vchg[i], lin, beta_tdc, cuts, eta);
      target += log_diff_exp(-HL, -HU) + Hent;    // {S(L)-S(U)}/S(A)
    }
  }
}

generated quantities {
  real sigma_out = use_frailty ? sigma_b[1] : 0;
}
