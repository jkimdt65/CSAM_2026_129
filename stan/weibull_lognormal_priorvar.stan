// weibull_lognormal_priorvar.stan
// Prior sensitivity for sigma_b (Section 5): the likelihood and the other
// priors of weibull_lognormal_frailty_ic.stan, with the prior on sigma_b
// selected by prior_id (1: half-N(0, 5^2), 2: half-N(0, 10^2),
// 3: exponential(0.5)).

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
  int<lower=1, upper=3> prior_id;   // 1: hN(0,5)  2: hN(0,10)  3: exp(0.5)
}
parameters {
  real<lower=0> sigma_b;
  real<lower=0> lambda;
  real<lower=0> gamma_shape;
  real beta;
  vector[n_clusters] z;
}
transformed parameters {
  vector[n_clusters] w = exp(sigma_b * z);
}
model {
  // Only the prior on sigma_b depends on prior_id; the likelihood below is the same for all.
  if (prior_id == 1) sigma_b ~ normal(0, 5);
  else if (prior_id == 2) sigma_b ~ normal(0, 10);
  else sigma_b ~ exponential(0.5);

  lambda      ~ normal(0, 5);
  gamma_shape ~ normal(0, 5);
  beta        ~ normal(0, 5);
  z           ~ std_normal();

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
  real hazard_ratio = exp(beta);
  real var_log_frailty = square(sigma_b);
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

