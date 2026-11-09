// sim_isp_tdc.stan
// Simulation of Section 4: I-spline baseline + time-varying covariate +
// delayed entry + (optional) log-normal frailty.
//
// Non-centred log-normal frailty:
//   log w_i = sigma_b * z_i,  z_i ~ N(0, 1)
// use_frailty = 0 fixes w = 1.
//
// Time-varying covariate: x2 switches from 0 to 1 at vchg (exogenous).
//   H(t) = w * exp(X*beta) * [ H0(min(t,v)) + (H0(t)-H0(v))_+ * exp(beta_tdc) ]
// Delayed entry: the cumulative hazard at entry is subtracted.

functions {
  // subject-level cumulative hazard with the time-varying covariate (basis computed in R)
  real cumhaz_tdc_isp(row_vector Ib_t, row_vector Ib_v,
                      real lin, real beta_tdc, vector gam, int after_v) {
    real Hpre  = dot_product(Ib_v, gam);
    real Hpost = after_v ? (dot_product(Ib_t, gam) - Hpre) : 0;
    real Ht    = after_v ? Hpre + Hpost * exp(beta_tdc)
                         : dot_product(Ib_t, gam);
    return exp(lin) * Ht;
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

  int<lower=2> nbasis;
  matrix[N, nbasis] IL;        // I-spline basis at L
  matrix[N, nbasis] IU;        // at U
  matrix[N, nbasis] Ivchg;     // at the covariate change point
  matrix[N, nbasis] Ient;      // at entry (zeros when entry = 0)
  array[N] int<lower=0, upper=1> L_after_v;
  array[N] int<lower=0, upper=1> U_after_v;
  array[N] int<lower=0, upper=1> ent_after_v;

  int<lower=1> n_clusters;
  array[N] int<lower=1> cluster;
  int<lower=0, upper=1> use_frailty;
}

parameters {
  vector[P] beta;
  real beta_tdc;
  real alpha;
  vector[nbasis - 1] xi_free;
  real<lower=0> smooth_sd;

  array[use_frailty ? 1 : 0] real<lower=0> sigma_b;
  vector[use_frailty ? n_clusters : 0] z;
}

transformed parameters {
  vector[nbasis] xi  = append_row(xi_free, -sum(xi_free));
  vector[nbasis] gam = exp(alpha + xi);   // positive by construction

  vector[use_frailty ? n_clusters : 0] log_w;
  if (use_frailty) log_w = sigma_b[1] * z;
}

model {
  beta      ~ normal(0, 5);
  beta_tdc  ~ normal(0, 5);
  alpha     ~ normal(0, 5);
  smooth_sd ~ normal(0, 1);
  for (k in 2:nbasis) xi[k] ~ normal(xi[k - 1], smooth_sd);

  if (use_frailty) {
    sigma_b[1] ~ normal(0, 1);
    z          ~ std_normal();
  }

  for (i in 1:N) {
    real wi  = use_frailty ? exp(log_w[cluster[i]]) : 1.0;
    real lin = P > 0 ? dot_product(X[i], beta) : 0.0;
    real HL  = wi * cumhaz_tdc_isp(IL[i], Ivchg[i], lin, beta_tdc, gam,
                                   L_after_v[i]);
    real Hent = entry[i] > 0
              ? wi * cumhaz_tdc_isp(Ient[i], Ivchg[i], lin, beta_tdc, gam,
                                    ent_after_v[i])
              : 0.0;

    if (cens_type[i] == 2) {
      target += -HL + Hent;                       // S(L)/S(A)
    } else {
      real HU = wi * cumhaz_tdc_isp(IU[i], Ivchg[i], lin, beta_tdc, gam,
                                    U_after_v[i]);
      target += log_diff_exp(-HL, -HU) + Hent;    // {S(L)-S(U)}/S(A)
    }
  }
}

generated quantities {
  real sigma_out = use_frailty ? sigma_b[1] : 0;
}
