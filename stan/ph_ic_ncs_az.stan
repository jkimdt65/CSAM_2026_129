// ph_ic_ncs_az.stan
// Royston-Parmar baseline, PH model, interval-censored data, no frailty.
//   log H0(t) = alpha + NS(log t) b       (Section 2.2, Royston and Parmar 2002)
// The basis is built in R without an intercept column and centred, so the
// intercept alpha is needed to estimate the level of log H0.
// rp_smooth = 1 puts an RW1 prior on b;
// rp_smooth = 0 gives b independent N(0, 5^2) priors, the usual unpenalised form.
// Prior scales passed as data.
//
// cens_type: 1 = exact, 2 = right, 3 = left, 4 = interval

data {
  int<lower=1> N;
  int<lower=0> P;
  matrix[N, P] X;
  array[N] int<lower=1, upper=4> cens_type;
  vector<lower=0>[N] entry;
  vector<lower=0>[N] L;
  int<lower=1> df;
  matrix[N, df] NS_L;
  matrix[N, df] NS_U;
  matrix[N, df] NS_ent;
  matrix[N, df] dNS_L;
  real<lower=0> beta_sd;
  real<lower=0> smooth_scale;
  int<lower=0, upper=1> rp_smooth;
}

parameters {
  vector[P] beta;
  real alpha;
  vector[df] b;
  real<lower=0> smooth_sd;      // used only when rp_smooth = 1
}

model {
  beta ~ normal(0, beta_sd);
  alpha ~ normal(0, 5);
  smooth_sd ~ normal(0, smooth_scale);
  if (rp_smooth == 1) {
    b[1] ~ normal(0, 5);
    for (k in 2:df) b[k] ~ normal(b[k - 1], smooth_sd);
  } else {
    b ~ normal(0, 5);
  }

  for (i in 1:N) {
    real eta = X[i] * beta;
    real ex = exp(eta);
    real logS_entry = entry[i] > 0 ? -ex * exp(alpha + dot_product(NS_ent[i], b)) : 0;
    if (cens_type[i] == 1) {
      real H_t = exp(alpha + dot_product(NS_L[i], b));
      real h0  = H_t * dot_product(dNS_L[i], b) / L[i];
      target += log(h0) + eta - ex * H_t - logS_entry;
    } else if (cens_type[i] == 2) {
      target += -ex * exp(alpha + dot_product(NS_L[i], b)) - logS_entry;
    } else if (cens_type[i] == 3) {
      target += log1m_exp(-ex * exp(alpha + dot_product(NS_U[i], b))) - logS_entry;
    } else {
      target += log_diff_exp(-ex * exp(alpha + dot_product(NS_L[i], b)),
                             -ex * exp(alpha + dot_product(NS_U[i], b))) - logS_entry;
    }
  }
}
