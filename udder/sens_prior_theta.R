# -----------------------------------------------------------------
# udder/sens_prior_theta.R
# Udder quarter infection data: sensitivity to the prior on the gamma
# frailty variance theta (Supplementary Table S2).
#
# The Weibull gamma frailty model of Table 2 is refitted with four priors
# on theta; the likelihood and all other priors are unchanged
# (weibull_gamma_priorvar.stan, prior selected by prior_id):
#   1: half-N(0, 5^2)          (the paper)
#   2: half-N(0, 10^2)
#   3: half-Cauchy(0, 1)
#   4: inverse gamma(0.01, 0.01), i.e. 1/theta ~ Gamma(0.01, 0.01),
#      the usual BUGS and JAGS choice for a frailty precision
# Setting 1 is refitted through the same program so that the four fits
# are comparable.
#
# Settings as in 02_fit_frailty.R: 4 chains x (1000 warmup + 4000
# sampling), seed 2026, adapt_delta 0.95. CmdStan 2.37.0. The four
# settings are fitted one after another.
#
# Run from the udder/ directory:
#   Rscript sens_prior_theta.R > sens_prior_theta_log.txt 2>&1
# Output: fit_prior_theta_{1..4}.rds, sens_prior_theta.md
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")

source("01_data_weibull.R")   # stan_data, d

OUTROOT <- "sens_theta_out"; dir.create(OUTROOT, showWarnings = FALSE)
lab <- c("half-N(0, 5^2)", "half-N(0, 10^2)", "half-Cauchy(0, 1)", "inv-gamma(0.01, 0.01)")

mod <- cmdstan_model(file.path("..", "stan", "weibull_gamma_priorvar.stan"))

run_one <- function(pid) {
  rds <- sprintf("fit_prior_theta_%d.rds", pid)
  cat(format(Sys.time(), "%H:%M:%S"), "prior", pid, ": start\n")
  od <- file.path(OUTROOT, paste0("p", pid)); dir.create(od, showWarnings = FALSE)
  t0 <- Sys.time()
  fit <- mod$sample(data = c(stan_data, list(prior_id = pid)),
                    chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 4000,
                    seed = 2026, adapt_delta = 0.95,
                    output_dir = od, refresh = 0, show_messages = FALSE)
  fit$save_object(rds)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(format(Sys.time(), "%H:%M:%S"), "prior", pid, ": done (", round(secs), "s)\n")
  list(fit = fit, pid = pid, secs = secs)
}

res <- lapply(1:4, function(pid) tryCatch(run_one(pid), error = function(e) {
  cat("prior", pid, "error:", conditionMessage(e), "\n"); NULL }))
stopifnot(!any(vapply(res, is.null, TRUE)))

f95 <- function(v) sprintf("%.3f (%.3f, %.3f)", median(v), quantile(v, .025), quantile(v, .975))
gd  <- function(f, p) as.numeric(f$draws(p, format = "draws_matrix"))

# Stan version, chain return codes (all 0) and number of draws (16000)
for (r in res) {
  m <- r$fit$metadata()
  cat("prior", r$pid, ": Stan", m$stan_version_major, m$stan_version_minor,
      m$stan_version_patch, "/ return codes", r$fit$return_codes(),
      "/ draws", posterior::ndraws(r$fit$draws("lp__")), "\n")
}

S <- do.call(rbind, lapply(res, function(r) {
  f <- r$fit
  data.frame(prior = lab[r$pid],
             theta = f95(gd(f, "theta")), beta = f95(gd(f, "beta")),
             tau = f95(gd(f, "kendall_tau")),
             ndiv = sum(f$diagnostic_summary()$num_divergent),
             rhat = sprintf("%.3f", max(f$summary(c("theta", "lambda", "gamma_shape", "beta"))$rhat)),
             secs = as.integer(round(r$secs)),
             stringsAsFactors = FALSE)
}))
cat("\n---- Summary ----\n"); print(S, row.names = FALSE)

th <- vapply(res, function(r) median(gd(r$fit, "theta")), 0)
be <- vapply(res, function(r) median(gd(r$fit, "beta")), 0)
cat(sprintf("\ntheta medians %.3f to %.3f (relative range %.1f%% of setting 1)\n",
            min(th), max(th), 100 * (max(th) - min(th)) / th[1]))
cat(sprintf("beta medians %.3f to %.3f (largest difference %.3f)\n",
            min(be), max(be), max(be) - min(be)))

L <- c("# Prior sensitivity for theta (gamma frailty), udder quarter data", "",
       "Weibull baseline, covariate multiparity. Posterior median (95% CrI).",
       "Same likelihood in all four rows; only the prior on theta differs.",
       "4 x (1000 + 4000), seed 2026, adapt_delta 0.95, CmdStan 2.37.0.", "",
       "| Prior on theta | theta | beta | Kendall tau | Divergences | max Rhat | Seconds |",
       "|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(S))) L <- c(L, with(S[i, ],
  sprintf("| %s | %s | %s | %s | %d | %s | %d |", prior, theta, beta, tau, ndiv, rhat, secs)))
L <- c(L, "",
  sprintf("theta medians %.3f to %.3f (relative range %.1f%% of the paper's setting).",
          min(th), max(th), 100 * (max(th) - min(th)) / th[1]),
  sprintf("beta medians %.3f to %.3f.", min(be), max(be)))
writeLines(L, "sens_prior_theta.md")
cat("\n-> sens_prior_theta.md written\n")
