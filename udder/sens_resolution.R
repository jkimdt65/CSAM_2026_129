# -----------------------------------------------------------------
# udder/sens_resolution.R
# Udder quarter infection data: sensitivity of the gamma frailty model
# to the resolution of the flexible baselines (Supplementary Table S1).
#
# The question: Table 4 shows Var(log w) moving by about 19% between the
# Weibull, piecewise exponential and I-spline baselines. If changing only
# the resolution within one construction moves it as much, the movement
# reflects flexibility in general rather than the construction.
#
# Design: piecewise exponential with K = 5, 8, 12 intervals; I-spline with
# 3, 5, 8 interior knots (6 fits). K = 8 and 5 knots are the settings of
# Table 4; they are refitted here through the same code path so that all
# six fits are comparable. Cutpoint and knot rules as in
# 03_data_baselines.R.
#
# Settings as in 04_fit_baselines.R: 4 chains x (1000 warmup + 4000
# sampling), seed 2026, adapt_delta 0.95. CmdStan 2.37.0.
# Three fits run at a time with 4 chains each (12 cores). Models are
# compiled in the main process first, and each fit writes to its own
# output directory.
#
# Run from the udder/ directory, after 04_fit_baselines.R (the between-
# construction variation is computed from its three spec-A fits):
#   Rscript sens_resolution.R > sens_resolution_log.txt 2>&1
# Output: fit_sens_*.rds, sens_resolution.rds, sens_resolution.md
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
library(splines2)
library(parallel)

source("03_data_baselines.R")   # d, stan_data_wb_A, cov_names_A, make_pwe_cuts()

N_CHAINS  <- 4
MC_CORES  <- 3                  # 3 fits x 4 chains = 12 cores
SEED      <- 2026
ADAPT_D   <- 0.95
OUTROOT   <- "sens_out"
dir.create(OUTROOT, showWarnings = FALSE)

# ===================================================================
# 1. Stan data for each setting
# ===================================================================
base_extra <- list(entry = rep(0, nrow(d)))

make_pwe_data <- function(K) {
  cuts <- make_pwe_cuts(L = d$L_stan, U = d$U_stan, K = K, method = "turnbull")
  K2 <- length(cuts) - 1
  ev <- d$U_stan[d$cens_type == 4]
  cnt <- table(cut(ev, breaks = cuts, include.lowest = TRUE))
  list(data = c(stan_data_wb_A, base_extra, list(K = K2, cuts = cuts)),
       K = K2, min_events = min(cnt))
}

make_isp_data <- function(nk, degree = 3) {
  pos <- c(d$L_stan[d$L_stan > 0], d$U_stan)
  bnd <- c(0, max(pos) * 1.05)
  ik  <- quantile(pos, probs = seq_len(nk) / (nk + 1))
  ib  <- iSpline(d$L_stan, knots = ik, degree = degree,
                 intercept = TRUE, Boundary.knots = bnd)
  IL <- predict(ib, d$L_stan); IU <- predict(ib, d$U_stan)
  stopifnot(all(abs(IL[d$L_stan == 0, ]) < 1e-12))           # H0(0) = 0
  stopifnot(all(IU[d$cens_type == 4, ] - IL[d$cens_type == 4, ] >= -1e-10))
  list(data = c(stan_data_wb_A, base_extra,
                list(nbasis = ncol(IL), IL = IL, IU = IU,
                     Ient = matrix(0, nrow(d), ncol(IL)))),
       nbasis = ncol(IL), min_events = NA)
}

cfg <- list()
for (K in c(5, 8, 12)) {
  m <- make_pwe_data(K)
  cfg[[length(cfg) + 1]] <- list(
    tag = sprintf("pwe_K%02d", m$K), family = "PWE", res = m$K,
    stan = "ph_ic_pwe_gamma_frailty.stan", data = m$data,
    note = sprintf("K = %d (min. %d events per interval)", m$K, m$min_events))
}
for (nk in c(3, 5, 8)) {
  m <- make_isp_data(nk)
  cfg[[length(cfg) + 1]] <- list(
    tag = sprintf("isp_k%d", nk), family = "I-spline", res = nk,
    stan = "ph_ic_ispline_gamma_frailty.stan", data = m$data,
    note = sprintf("%d interior knots (%d basis functions)", nk, m$nbasis))
}

cat("---- Settings ----\n")
for (c1 in cfg) cat(sprintf("%-10s %-9s %s\n", c1$tag, c1$family, c1$note))

# ===================================================================
# 2. Compile in the main process (concurrent compilation inside forked
#    workers would collide on temporary files)
# ===================================================================
cat("\n---- Compiling ----\n")
mods <- list()
for (f in unique(vapply(cfg, function(x) x$stan, ""))) {
  cat(f, "...\n"); mods[[f]] <- cmdstan_model(file.path("..", "stan", f))
}

# ===================================================================
# 3. Fits, three at a time. Each worker saves its fit and returns a
#    one-row summary only.
# ===================================================================
run_one <- function(k) {
  c1 <- cfg[[k]]
  od <- file.path(OUTROOT, c1$tag); dir.create(od, showWarnings = FALSE)
  t0 <- Sys.time()
  fit <- mods[[c1$stan]]$sample(
    data = c1$data, chains = N_CHAINS, parallel_chains = N_CHAINS,
    iter_warmup = 1000, iter_sampling = 4000,
    seed = SEED, adapt_delta = ADAPT_D,
    output_dir = od, refresh = 0, show_messages = FALSE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fit$save_object(paste0("fit_sens_", c1$tag, ".rds"))

  s  <- fit$summary(c("theta", "var_log_frailty", "beta"))
  dr <- fit$draws(c("alpha", "mean_log_w"), format = "draws_df")
  q  <- function(par, idx = 1) {
    dm <- fit$draws(par, format = "draws_matrix")[, idx]
    c(median(dm), quantile(dm, .025), quantile(dm, .975))
  }
  th <- q("theta"); vw <- q("var_log_frailty")
  bf <- q("beta", 1); bp6 <- q("beta", 7)   # front, parity 6
  data.frame(
    tag = c1$tag, family = c1$family, res = c1$res, note = c1$note,
    theta = th[1], theta_lo = th[2], theta_hi = th[3],
    varlogw = vw[1], varlogw_lo = vw[2], varlogw_hi = vw[3],
    front = bf[1], front_lo = bf[2], front_hi = bf[3],
    parity6 = bp6[1],
    ndiv = sum(fit$diagnostic_summary()$num_divergent),
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess = min(s$ess_bulk, na.rm = TRUE),
    ridge = cor(dr$alpha, dr$mean_log_w),
    secs = el, stringsAsFactors = FALSE)
}

cat("\n---- Fitting (", MC_CORES, " at a time) ----\n", sep = "")
t_all <- system.time(
  res <- mclapply(seq_along(cfg), run_one,
                  mc.cores = MC_CORES, mc.preschedule = FALSE)
)
bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), TRUE)
if (any(bad)) {
  cat("failed:", paste(vapply(cfg[bad], function(x) x$tag, ""), collapse = ", "), "\n")
  print(res[bad])
}
stopifnot(!any(bad))
R <- do.call(rbind, res)
cat("total elapsed:", round(t_all[["elapsed"]] / 60, 1), "min\n")

# ===================================================================
# 4. Within- vs between-construction variation in Var(log w)
# ===================================================================
# Relative variation = (max - min) / mean of the posterior medians.
# Between constructions: the three spec-A fits of Table 4 (04_fit_baselines.R).
rng <- function(v) 100 * (max(v) - min(v)) / mean(v)
vlw <- function(rds) median(as.numeric(readRDS(rds)$draws("var_log_frailty",
                                                         format = "draws_matrix")))
between_med <- c(Weibull = vlw("fit_wb_A.rds"), PWE = vlw("fit_pwe_A.rds"),
                 `I-spline` = vlw("fit_isp_A.rds"))
BETWEEN_PCT <- rng(between_med)
within <- tapply(R$varlogw, R$family, rng)
front_within <- tapply(R$front, R$family, rng)

cat("\n---- Relative variation in Var(log w) ----\n")
cat(sprintf("between constructions (Table 4)   : %.1f%%\n", BETWEEN_PCT))
for (f in names(within))
  cat(sprintf("within %-9s (resolution only) : %.1f%%\n", f, within[f]))
cat(sprintf("\nfront coefficient, within construction: PWE %.1f%%, I-spline %.1f%%\n",
            front_within["PWE"], front_within["I-spline"]))

cat("\n---- All fits ----\n")
print(R[, c("tag", "res", "theta", "varlogw", "front", "parity6",
            "ndiv", "max_rhat", "min_ess", "ridge", "secs")], digits = 4)

# ===================================================================
# 5. Summary (markdown)
# ===================================================================
f3 <- function(m, lo, hi) sprintf("%.3f (%.3f, %.3f)", m, lo, hi)
lines <- c(
  "# Baseline resolution (Supplementary Table S1)",
  "",
  "Gamma frailty, spec A (P = 8), 4,784 quarters, 1,196 cows.",
  "Posterior median (95% CrI). Seed and sampler settings as in Table 4.",
  "",
  "| Baseline | Resolution | theta | Var(log w) | front | parity6 | Divergences | max Rhat | min ESS | cor(alpha, mean log w) | Seconds |",
  "|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(R))) with(R[i, ], {
  lines <<- c(lines, sprintf("| %s | %s | %s | %s | %s | %.3f | %d | %.3f | %.0f | %+.3f | %.0f |",
    family, note, f3(theta, theta_lo, theta_hi),
    f3(varlogw, varlogw_lo, varlogw_hi), f3(front, front_lo, front_hi),
    parity6, ndiv, max_rhat, min_ess, ridge, secs)) })
lines <- c(lines, "",
  "## Relative variation in Var(log w)", "",
  sprintf("- between constructions (Weibull, PWE, I-spline; Table 4): %.1f%%", BETWEEN_PCT),
  sprintf("- within %s (resolution only): %.1f%%", names(within), within),
  "",
  sprintf("Front coefficient, within construction: PWE %.1f%%, I-spline %.1f%%.",
          front_within["PWE"], front_within["I-spline"]),
  "",
  sprintf("6 fits, 3 at a time x 4 chains, %.1f min in all.", t_all[["elapsed"]] / 60))
writeLines(lines, "sens_resolution.md")
saveRDS(list(R = R, between = between_med, BETWEEN_PCT = BETWEEN_PCT,
             within = within, front_within = front_within), "sens_resolution.rds")
cat("\n-> sens_resolution.md / .rds written\n")
