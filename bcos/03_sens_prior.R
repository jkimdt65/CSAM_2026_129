# -----------------------------------------------------------------
# bcos/03_sens_prior.R
# Breast cosmesis data: prior sensitivity (Supplementary Table S3).
#
# The priors of 02_fit_baselines.R are changed one at a time:
#   beta ~ N(0, s^2),            s = 2.5, 5 (paper), 10  (smoothing scale fixed at 1)
#   smoothing SD ~ half-N(0, c^2), c = 0.5, 1 (paper), 2  (beta SD fixed at 5)
# The Royston-Parmar fit has no random-walk smoothing, so only the beta
# prior is varied for it.
#
# Settings as in 02_fit_baselines.R: 4 chains x (1000 warmup + 1000
# sampling), adapt_delta 0.99, seed 2026. CmdStan 2.37.0.
#
# Run from the bcos/ directory:
#   Rscript 03_sens_prior.R > 03_sens_prior_log.txt 2>&1
# Output: bcos_sens_prior.md
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
library(parallel)

source("01_data.R")

# Royston-Parmar initial values, as in 02_fit_baselines.R
grid_t  <- exp(seq(log(min(d$U_stan[d$U_stan < BIG])), log(max(pos_rp)), length.out = 200))
NS_grid <- sweep(predict(ns_basis, log(grid_t)), 2, ctr, "-")
b0 <- as.numeric(qr.solve(NS_grid, log(grid_t) - mean(log(grid_t))))
init_rp <- function() list(alpha = -1.5, b = b0)

mods <- list(pwe = cmdstan_model(file.path("..", "stan", "ph_ic_pwe_az.stan")),
             isp = cmdstan_model(file.path("..", "stan", "ph_ic_ispline_az.stan")),
             rp  = cmdstan_model(file.path("..", "stan", "ph_ic_ncs_az.stan")))
base <- list(pwe = stan_data_pwe_bcos, isp = stan_data_isp_bcos,
             rp  = c(stan_data_rp_bcos, list(rp_smooth = 0L)))

grid <- rbind(
  expand.grid(model = c("pwe", "isp", "rp"), beta_sd = c(2.5, 5, 10), smooth_scale = 1,
              stringsAsFactors = FALSE),
  expand.grid(model = c("pwe", "isp"), beta_sd = 5, smooth_scale = c(0.5, 2),
              stringsAsFactors = FALSE))

OUTROOT <- "sens_bcos_out"; dir.create(OUTROOT, showWarnings = FALSE)

run_one <- function(k) {
  g <- grid[k, ]
  od <- file.path(OUTROOT, sprintf("run%02d", k)); dir.create(od, showWarnings = FALSE)
  fit <- mods[[g$model]]$sample(
    data = c(base[[g$model]], list(beta_sd = g$beta_sd, smooth_scale = g$smooth_scale)),
    chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000,
    seed = 2026, adapt_delta = 0.99, refresh = 0, show_messages = FALSE,
    output_dir = od,
    init = if (g$model == "rp") init_rp else NULL)
  v <- as.numeric(fit$draws("beta[1]", format = "draws_matrix"))
  data.frame(model = g$model, beta_sd = g$beta_sd, smooth_scale = g$smooth_scale,
             beta = sprintf("%.3f (%.3f, %.3f)", median(v), quantile(v, .025), quantile(v, .975)),
             med = median(v),
             ndiv = sum(fit$diagnostic_summary()$num_divergent),
             ok = all(fit$return_codes() == 0) && posterior::ndraws(fit$draws("lp__")) == 4000,
             stringsAsFactors = FALSE)
}

S <- do.call(rbind, mclapply(seq_len(nrow(grid)), run_one, mc.cores = 3, mc.preschedule = FALSE))
S <- S[order(S$model, S$smooth_scale, S$beta_sd), ]
cat("\n---- Summary ----\n"); print(S[, c("model", "beta_sd", "smooth_scale", "beta", "ndiv", "ok")], row.names = FALSE)
stopifnot(all(S$ok))

lab <- c(pwe = "Piecewise exponential", isp = "I-spline", rp = "Royston-Parmar")
L <- c("# Prior sensitivity, breast cosmesis data (revised models)", "",
       "Treatment coefficient, posterior median (95% CrI). adapt_delta 0.99, CmdStan 2.37.0, seed 2026.",
       "Paper settings: beta SD 5, smoothing scale 1; one prior changed at a time.", "",
       "| Baseline | beta SD | Smoothing scale | beta (RadChem) | Divergences |", "|---|---|---|---|---|")
for (i in seq_len(nrow(S))) L <- c(L, sprintf("| %s | %g | %s | %s | %d |",
  lab[S$model[i]], S$beta_sd[i],
  if (S$model[i] == "rp") "n/a" else format(S$smooth_scale[i]), S$beta[i], S$ndiv[i]))
rng <- tapply(S$med, S$model, function(v) max(v) - min(v))
L <- c(L, "", paste0("Largest change in the posterior median of beta: ",
  paste(sprintf("%s %.3f", lab[names(rng)], rng), collapse = ", "), "."))
writeLines(L, "bcos_sens_prior.md")
cat("\n-> bcos_sens_prior.md written\n")
