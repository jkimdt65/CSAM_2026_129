# -----------------------------------------------------------------
# udder/02_fit_frailty.R
# Udder quarter infection data, Weibull baseline: gamma and log-normal
# frailty models fitted to the same data with the same seed.
#   gamma frailty      : replication of Goethals et al. (2009), Table 2
#   gamma vs log-normal: Table 3 (Var(log w), Kendall's tau, LOO)
#
# Settings: 4 chains x (1000 warmup + 4000 sampling), seed 2026,
#           adapt_delta 0.95 (gamma) and 0.99 (log-normal). CmdStan 2.37.0.
#
# Run from the udder/ directory:
#   Rscript 02_fit_frailty.R > 02_fit_frailty_log.txt 2>&1
# Output: fit_gamma.rds, fit_lognormal.rds, loo_frailty.rds,
#         frailty_comparison.md
# The 95% intervals reported in Tables 2 and 3 are produced from the
# saved fits by 05_tables_95.R.
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
library(loo)

source("01_data_weibull.R")   # stan_data, d

N_CHAINS <- 4
ADAPT_D  <- 0.95     # gamma frailty; the log-normal fit uses 0.99 (below)
SEED     <- 2026

# ===================================================================
# 1. Gamma frailty (replication)
# ===================================================================
cat("\n========== gamma frailty ==========\n")
mod_g <- cmdstan_model("../stan/weibull_gamma_frailty_ic.stan")

t_g <- system.time(
  fit_g <- mod_g$sample(
    data            = stan_data,
    chains          = N_CHAINS,
    parallel_chains = N_CHAINS,
    iter_warmup     = 1000,
    iter_sampling   = 4000,
    seed            = SEED,
    adapt_delta     = ADAPT_D
  )
)
fit_g$cmdstan_diagnose()
fit_g$save_object("fit_gamma.rds")

# ===================================================================
# 2. Log-normal frailty
# ===================================================================
# adapt_delta is raised to 0.99 because, even in the non-centred form, a
# ridge remains between sigma_b and lambda: E[w] = exp(sigma_b^2 / 2)
# enters the likelihood only through its product with lambda.
#
# Initial values: the default initialisation on (-2, 2) of the
# unconstrained scale allows sigma_b up to e^2 = 7.4 and z up to 2, i.e.
# frailties up to exp(14.8). This does not produce NaN, since the
# likelihood is evaluated on the log scale, but it wastes the early
# warmup, so the chains start from moderate values.
cat("\n========== log-normal frailty ==========\n")
mod_l <- cmdstan_model("../stan/weibull_lognormal_frailty_ic.stan")

init_fun <- function() list(
  sigma_b     = 1.0,
  lambda      = 0.1,
  gamma_shape = 2.0,
  beta        = 0.5,
  z           = rnorm(stan_data$n_clusters, 0, 0.1)
)

t_l <- system.time(
  fit_l <- mod_l$sample(
    data            = stan_data,
    chains          = N_CHAINS,
    parallel_chains = N_CHAINS,
    iter_warmup     = 1000,
    iter_sampling   = 4000,
    seed            = SEED,
    adapt_delta     = 0.99,
    init            = init_fun
  )
)
fit_l$cmdstan_diagnose()
fit_l$save_object("fit_lognormal.rds")

cat("\nElapsed seconds: gamma", round(t_g[["elapsed"]], 1),
    "/ log-normal", round(t_l[["elapsed"]], 1), "\n")

# ===================================================================
# 3. Parameter summaries
# ===================================================================
sum_g <- fit_g$summary(c("theta", "lambda", "gamma_shape", "beta",
                         "kendall_tau", "hazard_ratio"))
sum_l <- fit_l$summary(c("sigma_b", "lambda", "lambda_mean1", "gamma_shape",
                         "beta", "var_log_frailty", "hazard_ratio"))
cat("\n--- gamma ---\n");      print(sum_g)
cat("\n--- log-normal ---\n"); print(sum_l)

# Comparison with the published maximum likelihood estimates, obtained by
# analytic integration of the gamma frailty (Goethals et al., 2009)
mle <- data.frame(
  parameter = c("theta", "lambda", "gamma_shape", "beta"),
  mle_est   = c(3.7967246, 0.1201593, 1.9672298, 0.8590531),
  mle_se    = sqrt(c(0.049281911, 0.0001982213, 0.0017306214, 0.020269244))
)
mle$mle_lo <- mle$mle_est - 1.96 * mle$mle_se
mle$mle_hi <- mle$mle_est + 1.96 * mle$mle_se
cat("\n--- gamma frailty: published MLE and posterior median ---\n")
print(merge(mle, as.data.frame(sum_g)[, c("variable", "median")],
            by.x = "parameter", by.y = "variable"), row.names = FALSE)

# ===================================================================
# 4. Common scale 1: variance of the log frailty
# ===================================================================
# gamma      : Var(log w) = trigamma(1/theta)   (closed form)
# log-normal : Var(log w) = sigma_b^2
# theta and sigma_b^2 are not on the same scale and are not compared directly.
draws_g <- fit_g$draws(format = "draws_df")
draws_l <- fit_l$draws(format = "draws_df")

vlog_g <- trigamma(1 / draws_g$theta)
vlog_l <- draws_l$var_log_frailty

qs <- function(v) c(median = median(v), q2.5 = quantile(v, .025),
                    q97.5 = quantile(v, .975))
cat("\nVar(log w)  gamma      :", round(qs(vlog_g), 3), "\n")
cat("Var(log w)  log-normal :", round(qs(vlog_l), 3), "\n")

# ===================================================================
# 5. Common scale 2: Kendall's tau
# ===================================================================
# For the gamma, tau = theta / (theta + 2); the log-normal has no closed
# form. The general expression (Hougaard, 2000) is
#   tau = 4 * int_0^inf s L(s) L''(s) ds - 1,
#   L(s) = E[exp(-sW)], L''(s) = E[W^2 exp(-sW)].
# For the log-normal the expectations are computed by Gauss-Hermite
# quadrature and the outer integral, after substituting s = exp(u), by
# adaptive quadrature.

# Gauss-Hermite nodes (probabilists' convention, weights summing to 1),
# by the Golub-Welsch eigendecomposition of the Jacobi matrix
gauss_hermite_prob <- function(n) {
  k <- 1:(n - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- sqrt(k)
  J[cbind(k + 1, k)] <- sqrt(k)
  e   <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = (e$vectors[1, ord])^2)
}
GH <- gauss_hermite_prob(120)
stopifnot(abs(sum(GH$w) - 1)            < 1e-10,   # weights sum to 1
          abs(sum(GH$w * GH$x^2) - 1)   < 1e-10,   # E[Z^2] = 1
          abs(sum(GH$w * GH$x^4) - 3)   < 1e-8)    # E[Z^4] = 3

# General tau from a Laplace transform and its second derivative.
# The range u in [-30, 30] corresponds to s in [1e-13, 1e13] and is stable
# to four decimal places for sigma_b up to 6.
kendall_tau_LT <- function(Lfun, L2fun, lo = -30, hi = 30) {
  g <- function(u) vapply(u, function(uu) {
    s <- exp(uu)
    s * s * Lfun(s) * L2fun(s)          # one factor s comes from ds = s du
  }, numeric(1))
  4 * integrate(g, lo, hi, subdivisions = 600L, rel.tol = 1e-10)$value - 1
}

# Check of the quadrature against the closed form for the gamma
tau_gamma_quad <- function(theta) {
  Lf  <- function(s) (1 + theta * s)^(-1 / theta)
  L2f <- function(s) (1 + theta) * (1 + theta * s)^(-1 / theta - 2)
  kendall_tau_LT(Lf, L2f)
}
cat("\n--- Quadrature check (gamma: quadrature vs closed form) ---\n")
for (th in c(0.5, 1, 3.797, 8)) {
  cat(sprintf("theta=%6.3f  quadrature=%.8f  closed form=%.8f\n",
              th, tau_gamma_quad(th), th / (th + 2)))
}

tau_lognormal <- function(sigma) {
  wv  <- exp(sigma * GH$x)
  Lf  <- function(s) sum(GH$w * exp(-s * wv))
  L2f <- function(s) sum(GH$w * wv^2 * exp(-s * wv))
  kendall_tau_LT(Lf, L2f)
}

# Applied to every fourth posterior draw (one quadrature per draw)
THIN     <- 4
sig_thin <- draws_l$sigma_b[seq(1, length(draws_l$sigma_b), by = THIN)]
cat("\nlog-normal tau: ", length(sig_thin), " draws ...\n", sep = "")
tau_l <- vapply(sig_thin, tau_lognormal, numeric(1))
tau_g <- draws_g$kendall_tau

cat("Kendall tau  gamma      :", round(qs(tau_g), 4), "\n")
cat("Kendall tau  log-normal :", round(qs(tau_l), 4), "\n")

# ===================================================================
# 6. Approximate leave-one-out comparison
# ===================================================================
# log_lik is the log-likelihood conditional on the sampled frailties, so
# (a) observation-level LOO leaves out one quarter while the frailty is
#     still informed by it, which is optimistic and produces high Pareto k;
# (b) summing log_lik within cows does not give a true leave-one-cow-out
#     either, since the frailty of the omitted cow is still learned from
#     its data.
# The same approximation is applied to both models, so the comparison is
# indicative only (Section 5 of the paper).
ll_g <- fit_g$draws("log_lik", format = "matrix")
ll_l <- fit_l$draws("log_lik", format = "matrix")

chain_id <- rep(1:N_CHAINS, each = nrow(ll_g) / N_CHAINS)

r_eff_g <- relative_eff(exp(ll_g), chain_id = chain_id)
r_eff_l <- relative_eff(exp(ll_l), chain_id = chain_id)

loo_g <- loo(ll_g, r_eff = r_eff_g)
loo_l <- loo(ll_l, r_eff = r_eff_l)

cat("\n--- observation-level LOO: gamma ---\n");      print(loo_g)
cat("\n--- observation-level LOO: log-normal ---\n"); print(loo_l)

cat("\nobservations with Pareto k > 0.7: gamma",
    sum(pareto_k_values(loo_g) > 0.7), " / log-normal",
    sum(pareto_k_values(loo_l) > 0.7), "\n")

cmp_obs <- loo_compare(list(gamma = loo_g, lognormal = loo_l))
cat("\n--- observation-level comparison ---\n"); print(cmp_obs)

# Summed within cows (from the log_lik matrix alone)
agg <- function(ll, cl) {
  out <- matrix(0, nrow(ll), max(cl))
  for (k in seq_len(max(cl))) out[, k] <- rowSums(ll[, cl == k, drop = FALSE])
  out
}
llc_g <- agg(ll_g, d$cluster_id)
llc_l <- agg(ll_l, d$cluster_id)

loo_cg <- loo(llc_g, r_eff = relative_eff(exp(llc_g), chain_id = chain_id))
loo_cl <- loo(llc_l, r_eff = relative_eff(exp(llc_l), chain_id = chain_id))
cmp_clu <- loo_compare(list(gamma = loo_cg, lognormal = loo_cl))
cat("\n--- cow-level comparison (approximate) ---\n"); print(cmp_clu)

saveRDS(list(loo_g = loo_g, loo_l = loo_l, cmp_obs = cmp_obs, cmp_clu = cmp_clu,
             tau_l = tau_l, seconds = c(gamma = t_g[["elapsed"]],
                                        lognormal = t_l[["elapsed"]])),
        "loo_frailty.rds")

# ===================================================================
# 7. Summary table (markdown)
# ===================================================================
f3 <- function(v) sprintf("%.3f (%.3f, %.3f)",
                          median(v), quantile(v, .025), quantile(v, .975))
gq <- function(fit, par) f3(fit$draws(par, format = "draws_matrix")[, 1])

lines <- c(
  "# Gamma vs log-normal frailty (Weibull baseline, udder quarter data)",
  "",
  "Posterior median (95% CrI). 4,784 quarters, 1,196 cows.",
  "",
  "| | Gamma frailty | Log-normal frailty |",
  "|---|---|---|",
  paste0("| Weibull shape | ", gq(fit_g, "gamma_shape"), " | ",
         gq(fit_l, "gamma_shape"), " |"),
  paste0("| beta (multiparous) | ", gq(fit_g, "beta"), " | ",
         gq(fit_l, "beta"), " |"),
  paste0("| Hazard ratio | ", gq(fit_g, "hazard_ratio"), " | ",
         gq(fit_l, "hazard_ratio"), " |"),
  paste0("| Weibull scale (frailty mean 1) | ", gq(fit_g, "lambda"), " | ",
         gq(fit_l, "lambda_mean1"), " |"),
  paste0("| Frailty parameter | theta = ", gq(fit_g, "theta"),
         " | sigma_b = ", gq(fit_l, "sigma_b"), " |"),
  paste0("| Var(log w) | ", f3(vlog_g), " | ", f3(vlog_l), " |"),
  paste0("| Kendall tau | ", f3(tau_g), " | ", f3(tau_l), " |"),
  paste0("| Divergent transitions | ", sum(fit_g$diagnostic_summary()$num_divergent),
         " | ", sum(fit_l$diagnostic_summary()$num_divergent), " |"),
  paste0("| Elapsed seconds | ", round(t_g[["elapsed"]], 0), " | ",
         round(t_l[["elapsed"]], 0), " |"),
  "",
  "## LOO comparison",
  "",
  "```", paste(capture.output(print(cmp_obs)), collapse = "\n"), "```",
  "",
  "```", paste(capture.output(print(cmp_clu)), collapse = "\n"), "```",
  "",
  paste0("Pareto k > 0.7: gamma ", sum(pareto_k_values(loo_g) > 0.7),
         ", log-normal ", sum(pareto_k_values(loo_l) > 0.7),
         " (of ", ncol(ll_g), " observations)"),
  "",
  "log_lik is conditional on the sampled frailties, so the comparison is",
  "not an exact leave-one-cow-out and is indicative only."
)
writeLines(lines, "frailty_comparison.md")
cat("\n-> frailty_comparison.md written\n")
