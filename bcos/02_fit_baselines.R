# -----------------------------------------------------------------
# bcos/02_fit_baselines.R
# Breast cosmesis data: the treatment effect under the three baseline
# constructions of Section 2.2 (Section 5 of the paper).
#
#   piecewise exponential, I-spline : log coefficients = alpha + xi,
#                                     sum(xi) = 0, random walk on xi
#   Royston-Parmar                   : log H0 = alpha + NS b, with
#                                     independent N(0, 5^2) priors on b
#
# Settings: 4 chains x (1000 warmup + 1000 sampling), adapt_delta 0.99,
#           seed 2026, prior SD 5 for beta and scale 1 for the
#           smoothing SD. CmdStan 2.37.0.
#
# Run from the bcos/ directory:
#   Rscript 02_fit_baselines.R > 02_fit_baselines_log.txt 2>&1
# Output: fit_bcos_{pwe,isp,rp}.rds, bcos_baselines.md
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
library(parallel)

source("01_data.R")

N_CHAINS <- 4
SEED     <- 2026
ADAPT_D  <- 0.99
OUTROOT  <- "bcos_out"
dir.create(OUTROOT, showWarnings = FALSE)

# Royston-Parmar initial values. Random initial values for b can make
# log H0 decrease over part of the range, so that H0(L) > H0(U) for some
# observations and the log-likelihood is undefined. The chains therefore
# start from the b0 that reproduces, by least squares on the NS basis, a
# straight line in log t (an exponential baseline), which is monotone.
grid_t  <- exp(seq(log(min(d$U_stan[d$U_stan < BIG])), log(max(pos_rp)), length.out = 200))
NS_grid <- sweep(predict(ns_basis, log(grid_t)), 2, ctr, "-")
b0 <- as.numeric(qr.solve(NS_grid, log(grid_t) - mean(log(grid_t))))
stopifnot(all(diff(as.numeric(NS_grid %*% b0)) > 0))
init_rp <- function() list(alpha = -1.5, b = b0)

pri <- list(beta_sd = 5, smooth_scale = 1)

specs <- list(
  pwe = list(stan = "ph_ic_pwe_az.stan",     data = c(stan_data_pwe_bcos, pri)),
  isp = list(stan = "ph_ic_ispline_az.stan", data = c(stan_data_isp_bcos, pri)),
  rp  = list(stan = "ph_ic_ncs_az.stan",     data = c(stan_data_rp_bcos, pri, list(rp_smooth = 0L)),
             init = init_rp)
)

cat("\n---- Compiling ----\n")
mods <- list()
for (f in unique(vapply(specs, `[[`, "", "stan"))) {
  cat(f, "...\n")
  mods[[f]] <- cmdstan_model(file.path("..", "stan", f))
}

run_one <- function(nm) {
  od <- file.path(OUTROOT, nm); dir.create(od, showWarnings = FALSE)
  t0 <- Sys.time()
  fit <- mods[[specs[[nm]]$stan]]$sample(
    data = specs[[nm]]$data, chains = N_CHAINS, parallel_chains = N_CHAINS,
    iter_warmup = 1000, iter_sampling = 1000,
    seed = SEED, adapt_delta = ADAPT_D, init = specs[[nm]]$init,
    output_dir = od, refresh = 0, show_messages = FALSE)
  fit$save_object(paste0("fit_bcos_", nm, ".rds"))
  list(fit = fit, name = nm,
       elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

cat("\n---- Fitting ----\n")
res <- mclapply(names(specs), run_one, mc.cores = 3, mc.preschedule = FALSE)
names(res) <- vapply(res, function(x) x$name, "")
res <- res[names(specs)]

# Stan version, chain return codes (all 0) and number of draws (4000)
for (nm in names(res)) {
  fit <- res[[nm]]$fit; m <- fit$metadata()
  cat(nm, ": Stan", m$stan_version_major, m$stan_version_minor,
      m$stan_version_patch, "/ return codes", fit$return_codes(),
      "/ draws", posterior::ndraws(fit$draws("lp__")), "\n")
}

# ---- Summary ----
f3 <- function(v) sprintf("%.3f (%.3f, %.3f)", median(v), quantile(v, .025), quantile(v, .975))
gd <- function(fit, p) as.numeric(fit$draws(p, format = "draws_matrix")[, 1])
ndiv <- function(fit) sum(fit$diagnostic_summary()$num_divergent)
mrhat <- function(fit) max(fit$summary()$rhat, na.rm = TRUE)

# Royston-Parmar monotonicity: alpha only shifts the level, so it is
# judged from b alone, as the share of draws with log H0 increasing on the grid
mono <- function(fit) {
  b <- as.matrix(fit$draws("b", format = "draws_matrix"))
  mean(apply(b %*% t(NS_grid), 1, function(v) all(diff(v) > 0)))
}

lab <- c(pwe = "Piecewise exponential", isp = "I-spline", rp = "Royston-Parmar")
rows <- lapply(names(res), function(nm) {
  fit <- res[[nm]]$fit
  data.frame(model = lab[nm], beta = f3(gd(fit, "beta[1]")),
             ndiv = ndiv(fit), rhat = sprintf("%.4f", mrhat(fit)),
             mono = if (nm == "rp") sprintf("%.1f%%", 100 * mono(fit)) else "",
             secs = round(res[[nm]]$elapsed), stringsAsFactors = FALSE)
})
S <- do.call(rbind, rows)
cat("\n---- Summary ----\n"); print(S, row.names = FALSE)
cat("\nRoyston-Parmar alpha:", f3(gd(res$rp$fit, "alpha")), "\n")

L <- c("# Breast cosmesis data: treatment effect under three baselines", "",
       "Posterior median (95% CrI). 4 x 1000/1000, adapt_delta 0.99, seed 2026, CmdStan 2.37.0.", "",
       "| Baseline | beta (RadChem) | Divergences | max Rhat | RP monotone | Seconds |",
       "|---|---|---|---|---|---|")
for (i in seq_len(nrow(S))) {
  L <- c(L, sprintf("| %s | %s | %d | %s | %s | %d |", S$model[i], S$beta[i],
                    S$ndiv[i], S$rhat[i], S$mono[i], S$secs[i]))
}
L <- c(L, "", paste("Royston-Parmar alpha:", f3(gd(res$rp$fit, "alpha"))))
writeLines(L, "bcos_baselines.md")
cat("\n-> bcos_baselines.md written\n")
