// -----------------------------------------------------------------
// weibull_lognormal_frailty_ic.stan
//
// Log-normal frailty counterpart to weibull_gamma_frailty_ic.stan
//   H(t | x, w) = w * lambda * t^gamma * exp(x*beta)
//   w_i = exp(b_i),  b_i = sigma_b * z_i,  z_i ~ N(0, 1)
//
// This file is deliberately kept LINE-FOR-LINE parallel to the gamma
// version: the functions block, the data block, the likelihood loop and
// the generated-quantities loop are byte-identical. Everything that
// differs is confined to `parameters`, `transformed parameters` and four
// lines of `model`. That is the point of section 4.5 -- lifting the gamma
// restriction is a local edit, not a redesign, because the frailty enters
// the likelihood only through the scalar multiplier w[cluster[i]].
//
// The closed-form (inclusion-exclusion) implementation cannot make this
// edit at all: its analytic marginalization exists only for the gamma.
//
// NORMALIZATION -- read before comparing lambda across the two models.
//   gamma version:      E[w] = 1,       median(w) != 1
//   this version:       E[b] = 0  =>    median(w) = 1, E[w] = exp(sigma_b^2/2)
// So lambda is on a different footing in the two fits and the raw values
// are NOT comparable. `lambda_mean1` in generated quantities rescales this
// model's lambda onto the mean-one-frailty footing used by the gamma model.
// The E[b] = 0 convention is chosen to match coxme / frailtypack, so the
// fitted sigma_b can be cross-checked against those packages directly.
//
// Time scale: t is in QUARTER-YEARS (original days / 91.31), as in the
// gamma version, so the two fits are directly comparable.
//
// censoring coding (cens_type):
//   2 = right censored   -> quarter infection-free through time L
//   4 = interval censored -> infection occurred in (L, U]
// -----------------------------------------------------------------

functions {
  // Weibull cumulative baseline hazard, safe at t = 0.
  // (Identical to the gamma version. Retained verbatim: the 68 quarters
  // with L = 0 are left-censored, H0(0) = 0 makes log_diff_exp(0, -HU)
  // reduce to log{1 - S(U)}, and the hard 0 short-circuits the NaN
  // gradient that pow(0, gamma) would otherwise produce in reverse-mode
  // autodiff.)
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
  real<lower=0> sigma_b;                 // sd of the LOG frailty
                                         //   (replaces theta)
  real<lower=0> lambda;                  // Weibull scale
  real<lower=0> gamma_shape;             // Weibull shape
  real beta;                             // log hazard ratio for x
  vector[n_clusters] z;                  // standardized frailties
                                         //   (replaces w, now unconstrained)
}

transformed parameters {
  // Non-centered parametrization. Two reasons, both worth stating in the
  // paper:
  //   (1) Geometry. Each cow contributes only 4 quarters, most of them
  //       right-censored, so the data are weak at the cluster level --
  //       exactly the regime where the centered parametrization produces
  //       Neal's funnel and divergent transitions. Non-centering makes the
  //       posterior for z approximately independent of sigma_b.
  //   (2) Support. w = exp(.) is positive by construction, so no
  //       <lower=0> declaration and no boundary. This also removes the
  //       `gamma_lpdf: Random variable is 0` initialization warnings that
  //       the gamma version emits (1-6 per chain), since nothing here can
  //       be initialized at zero.
  vector[n_clusters] w = exp(sigma_b * z);
}

model {
  // ---- priors ----
  // Matching the gamma version's weakly-informative choices in spirit.
  // Scale check: the replication gives theta = 3.80, and for a gamma
  // frailty Var(log w) = trigamma(1/theta) = 15.59, i.e. an equivalent
  // log-scale sd of about 3.95. A half-normal(0, 5) therefore places the
  // plausible region near 0.8 prior sd -- weakly informative, not
  // restrictive, but worth a sensitivity check (see the .md notes).
  sigma_b     ~ normal(0, 5);            // half-normal via <lower=0>
  lambda      ~ normal(0, 5);
  gamma_shape ~ normal(0, 5);
  beta        ~ normal(0, 5);

  z ~ std_normal();                      // <- the only structural change

  // ---- likelihood ----
  // Byte-identical to the gamma version below this line.
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
  real hazard_ratio = exp(beta);

  // ---- common footing with the gamma model ----
  // Kendall's tau has NO closed form for a log-normal frailty, so it is
  // computed by quadrature in R (compare_frailty_loo.R) rather than here.
  // What IS available in closed form is the variance of the frailty on the
  // log-hazard scale, which is the natural common currency:
  //   log-normal : Var(log w) = sigma_b^2
  //   gamma      : Var(log w) = trigamma(1/theta)
  // Note these two footings do NOT agree -- matching Var(log w) and
  // matching Kendall's tau give different sigma_b (3.95 vs 3.02 at
  // theta = 3.80). Report both and say so; the disagreement is
  // informative, not an error.
  real var_log_frailty = square(sigma_b);

  // lambda rescaled onto the mean-one-frailty convention of the gamma
  // model, so that Table 2 can put the two scale parameters side by side.
  real lambda_mean1 = lambda * exp(square(sigma_b) / 2);

  // ---- pointwise log-likelihood for LOO ----
  // Identical construction to the gamma version, so the two log_lik
  // matrices are directly comparable observation by observation.
  // CAVEAT (carried into the R script): this is the CONDITIONAL
  // log-likelihood given the sampled frailty. It supports
  // leave-one-observation-out, not leave-one-cow-out.
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
