# -----------------------------------------------------------------
# udder/sens_prior_sigma.R
# Udder quarter infection data: sensitivity of the log-normal frailty fit
# to the prior on sigma_b (Section 5).
#
# The question: Section 5 reports that the gamma frailty is favoured over
# the log-normal on the approximate LOO scale. If that depended on the
# prior for sigma_b it would be weak evidence. Three priors are compared:
#   A. sigma_b ~ half-N(0, 5^2)   the paper; the fit of 02_fit_frailty.R
#                                  (fit_lognormal.rds) is used as it is
#   B. sigma_b ~ half-N(0, 10^2)
#   C. sigma_b ~ exponential(0.5) (mean 2)
# B and C use weibull_lognormal_priorvar.stan, which has the likelihood
# and the other priors of weibull_lognormal_frailty_ic.stan and selects
# the prior on sigma_b through prior_id.
#
# Settings as in 02_fit_frailty.R for the log-normal model: 4 chains x
# (1000 warmup + 4000 sampling), seed 2026, adapt_delta 0.99, the same
# initial values. CmdStan 2.37.0. B and C run at the same time (8 cores).
#
# Run from the udder/ directory, after 02_fit_frailty.R:
#   Rscript sens_prior_sigma.R > sens_prior_sigma_log.txt 2>&1
# Output: fit_prior_B.rds, fit_prior_C.rds, sens_prior_sigma.md
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
library(loo)
library(parallel)

source("01_data_weibull.R")   # stan_data, d

SEED <- 2026; N_CHAINS <- 4
stopifnot(file.exists("fit_lognormal.rds"), file.exists("loo_frailty.rds"))

mod <- cmdstan_model(file.path("..", "stan", "weibull_lognormal_priorvar.stan"))

# ===================================================================
# 1. Fits B and C
# ===================================================================
init_fun <- function() list(sigma_b = 1.0, lambda = 0.1, gamma_shape = 2.0,
                            beta = 0.5, z = rnorm(stan_data$n_clusters, 0, 0.1))

run <- function(pid) {
  tag <- c("A", "B", "C")[pid]
  od  <- paste0("prior_out_", tag); dir.create(od, showWarnings = FALSE)
  fit <- mod$sample(data = c(stan_data, list(prior_id = pid)),
                    chains = N_CHAINS, parallel_chains = N_CHAINS,
                    iter_warmup = 1000, iter_sampling = 4000,
                    seed = SEED, adapt_delta = 0.99, init = init_fun,
                    output_dir = od, refresh = 0, show_messages = FALSE)
  fit$save_object(paste0("fit_prior_", tag, ".rds"))
  fit
}

cat("---- Fitting B (half-N(0, 10^2)) and C (exponential(0.5)) ----\n")
tt <- system.time(fits <- mclapply(2:3, run, mc.cores = 2, mc.preschedule = FALSE))
cat("done:", round(tt[["elapsed"]]), "s\n")
names(fits) <- c("B", "C")
fits$A <- readRDS("fit_lognormal.rds")
fits <- fits[c("A", "B", "C")]
lab  <- c(A = "half-N(0, 5^2)", B = "half-N(0, 10^2)", C = "exponential(0.5)")

# ===================================================================
# 2. Summary and LOO
# ===================================================================
f3 <- function(v) sprintf("%.3f (%.3f, %.3f)",
                          median(v), quantile(v, .025), quantile(v, .975))
gd <- function(f, p) as.numeric(f$draws(p, format = "draws_matrix")[, 1])

chain_id <- rep(1:N_CHAINS, each = 4000)
loos <- lapply(fits, function(f) {
  ll <- f$draws("log_lik", format = "matrix")
  loo(ll, r_eff = relative_eff(exp(ll), chain_id = chain_id))
})
cmp <- loo_compare(loos)

S <- do.call(rbind, lapply(names(fits), function(k) data.frame(
  prior   = lab[k],
  sigma_b = f3(gd(fits[[k]], "sigma_b")),
  beta    = f3(gd(fits[[k]], "beta")),
  elpd    = sprintf("%.1f", loos[[k]]$estimates["elpd_loo", "Estimate"]),
  ndiv    = sum(fits[[k]]$diagnostic_summary()$num_divergent),
  rhat    = sprintf("%.3f", max(fits[[k]]$summary(c("sigma_b", "lambda", "gamma_shape", "beta"))$rhat)),
  stringsAsFactors = FALSE)))
cat("\n---- Summary ----\n"); print(S, row.names = FALSE)

sig <- vapply(fits, function(f) median(gd(f, "sigma_b")), 0)
bet <- vapply(fits, function(f) median(gd(f, "beta")), 0)
el  <- vapply(loos, function(l) l$estimates["elpd_loo", "Estimate"], 0)
GAMMA_ELPD <- readRDS("loo_frailty.rds")$loo_g$estimates["elpd_loo", "Estimate"]
gap <- GAMMA_ELPD - max(el)

cat(sprintf("\nsigma_b medians %.3f to %.3f (relative range %.1f%%)\n",
            min(sig), max(sig), 100 * (max(sig) - min(sig)) / mean(sig)))
cat(sprintf("beta medians within %.3f of setting A\n", max(abs(bet - bet["A"]))))
cat(sprintf("elpd from %.1f to %.1f; gamma frailty %.1f; smallest gap %.1f\n",
            min(el), max(el), GAMMA_ELPD, gap))
cat("\n---- LOO comparison among the three priors ----\n"); print(cmp)

# ===================================================================
# 3. Summary (markdown)
# ===================================================================
L <- c("# Prior sensitivity for sigma_b (log-normal frailty), udder quarter data", "",
  "Weibull baseline, covariate multiparity. Posterior median (95% CrI).",
  "Same likelihood in all three rows; only the prior on sigma_b differs.", "",
  "| Prior on sigma_b | sigma_b | beta | elpd_loo | Divergences | max Rhat |",
  "|---|---|---|---|---|---|")
for (i in seq_len(nrow(S))) {
  L <- c(L, sprintf("| %s | %s | %s | %s | %d | %s |",
                    S$prior[i], S$sigma_b[i], S$beta[i], S$elpd[i], S$ndiv[i], S$rhat[i]))
}
L <- c(L, "",
  sprintf("sigma_b medians %.3f to %.3f (relative range %.1f%%); beta within %.3f of setting A.",
          min(sig), max(sig), 100 * (max(sig) - min(sig)) / mean(sig), max(abs(bet - bet["A"]))),
  sprintf("elpd from %.1f to %.1f. The gamma frailty has %.1f, so the smallest gap is %.1f.",
          min(el), max(el), GAMMA_ELPD, gap),
  "", "```", paste(capture.output(print(cmp)), collapse = "\n"), "```")
writeLines(L, "sens_prior_sigma.md")
cat("\n-> sens_prior_sigma.md written\n")
