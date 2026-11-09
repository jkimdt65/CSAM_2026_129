// weibull_gamma_priorvar.stan
// Prior sensitivity for the gamma frailty variance theta
// (Supplementary Table S2).
// The likelihood is identical to weibull_gamma_frailty_ic.stan (Table 2);
// only the prior on theta changes, selected by prior_id:
//   1: theta ~ half-N(0, 5^2)        (as in the manuscript)
//   2: theta ~ half-N(0, 10^2)
//   3: theta ~ half-Cauchy(0, 1)
//   4: theta ~ inv-gamma(0.01, 0.01), i.e. 1/theta ~ Gamma(0.01, 0.01),
//      the conventional BUGS/JAGS choice for a frailty precision
functions {
  real weib_cumhaz(real t, real lambda, real gamma_shape) {
    if (t <= 0) return 0;
    return lambda * pow(t, gamma_shape);
  }
}
data {
  int<lower=1> N;
  int<lower=1> n_clusters;
  array[N] int<lower=1> cluster;
  array[N] int<lower=2, upper=4> cens_type;
  vector<lower=0>[N] L;
  vector<lower=0>[N] U;
  vector[N] x;
  int<lower=1, upper=4> prior_id;
}
parameters {
  real<lower=0> theta;
  real<lower=0> lambda;
  real<lower=0> gamma_shape;
  real beta;
  vector<lower=0>[n_clusters] w;
}
model {
  if (prior_id == 1)      theta ~ normal(0, 5);
  else if (prior_id == 2) theta ~ normal(0, 10);
  else if (prior_id == 3) theta ~ cauchy(0, 1);
  else                    theta ~ inv_gamma(0.01, 0.01);

  lambda      ~ normal(0, 5);
  gamma_shape ~ normal(0, 5);
  beta        ~ normal(0, 5);
  w ~ gamma(1 / theta, 1 / theta);

  for (i in 1:N) {
    real ex = exp(x[i] * beta);
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
  real kendall_tau = theta / (theta + 2);
  real hazard_ratio = exp(beta);
}
