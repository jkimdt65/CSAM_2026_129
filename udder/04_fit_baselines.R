# -----------------------------------------------------------------
# udder/04_fit_baselines.R
# Udder quarter infection data: the gamma frailty model under three
# baseline constructions (Table 4), and the within- vs between-cluster
# comparison of posterior precision (Section 5).
#   1) Weibull                 + gamma frailty, spec A (P = 8)
#   2) piecewise exp. (K = 8)  + gamma frailty, spec A
#   3) I-spline (5 knots)      + gamma frailty, spec A
#   4) Weibull                 + gamma frailty, spec B (P = 3)
#
# Settings: 4 chains x (1000 warmup + 4000 sampling), seed 2026,
#           adapt_delta 0.95. CmdStan 2.37.0.
#
# Run from the udder/ directory:
#   Rscript 04_fit_baselines.R > 04_fit_baselines_log.txt 2>&1
# Output: fit_wb_A.rds, fit_pwe_A.rds, fit_isp_A.rds, fit_wb_B.rds,
#         baselines.md
# The 95% intervals of Table 4 are produced by 05_tables_95.R.
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")

source("03_data_baselines.R")

N_CHAINS <- 4
SEED     <- 2026
ADAPT_D  <- 0.95

# ===================================================================
# 1. Fits
# ===================================================================
fit_one <- function(stan_file, data, tag, adapt_delta = ADAPT_D) {
  rds <- paste0("fit_", tag, ".rds")
  cat("\n==========", tag, "==========\n")
  mod <- cmdstan_model(file.path("..", "stan", stan_file))
  tt <- system.time(
    fit <- mod$sample(
      data = data, chains = N_CHAINS, parallel_chains = N_CHAINS,
      iter_warmup = 1000, iter_sampling = 4000,
      seed = SEED, adapt_delta = adapt_delta
    )
  )
  fit$cmdstan_diagnose()
  fit$save_object(rds)
  cat(tag, "elapsed:", round(tt[["elapsed"]], 1), "seconds\n")
  attr(fit, "elapsed") <- tt[["elapsed"]]
  fit
}

fit_wb_A  <- fit_one("weibull_gamma_frailty_X.stan",      stan_data_wb_A,  "wb_A")
fit_pwe_A <- fit_one("ph_ic_pwe_gamma_frailty.stan",      stan_data_pwe_A, "pwe_A")
fit_wb_B  <- fit_one("weibull_gamma_frailty_X.stan",      stan_data_wb_B,  "wb_B")
fit_isp_A <- fit_one("ph_ic_ispline_gamma_frailty.stan",  stan_data_isp_A, "isp_A")

fits <- list(Weibull = fit_wb_A, `Piecewise exponential` = fit_pwe_A,
             `I-spline` = fit_isp_A)

# ===================================================================
# 2. Separation of baseline level and frailty
# ===================================================================
# A posterior correlation between alpha and mean_log_w close to -1 would
# mean that the level of the baseline and the frailties are not separated.
ridge <- function(fit, tag) {
  dr <- fit$draws(c("alpha", "mean_log_w"), format = "draws_df")
  r <- cor(dr$alpha, dr$mean_log_w)
  cat(sprintf("%-22s cor(alpha, mean_log_w) = %+.3f%s\n", tag, r,
              if (abs(r) > 0.9) "   (weak separation)" else ""))
  r
}
cat("\n---- Level and frailty separation ----\n")
rho <- c(`Piecewise exponential` = ridge(fit_pwe_A, "Piecewise exponential"),
         `I-spline` = ridge(fit_isp_A, "I-spline"))

# ===================================================================
# 3. Within- vs between-cluster precision (spec B)
# ===================================================================
sB <- fit_wb_B$summary("beta")
sB$name  <- cov_names_B
sB$scale <- ifelse(sB$name %in% c("front", "left"), "within", "between")
cat("\n---- Within- vs between-cluster precision (spec B) ----\n")
print(sB[, c("name", "scale", "median", "sd", "q5", "q95", "ess_bulk")])

w_sd <- mean(sB$sd[sB$scale == "within"])
b_sd <- mean(sB$sd[sB$scale == "between"])
cat("\nmean posterior SD: within", round(w_sd, 4), "/ between", round(b_sd, 4),
    "/ ratio", round(b_sd / w_sd, 2), "\n")
cat("  for reference: category balance alone would give 1.06, clustering alone 2.00\n")

# ===================================================================
# 4. Summary table (markdown)
# ===================================================================
f3 <- function(v) sprintf("%.3f (%.3f, %.3f)",
                          median(v), quantile(v, .025), quantile(v, .975))
par_row <- function(fit, par, idx = NULL) {
  dm <- fit$draws(par, format = "draws_matrix")
  f3(if (is.null(idx)) dm[, 1] else dm[, idx])
}
ndiv <- function(fit) sum(fit$diagnostic_summary()$num_divergent)
elap <- function(fit) as.character(round(attr(fit, "elapsed")))

lines <- c(
  "# Baseline construction and frailty (gamma frailty, spec A)",
  "",
  "Posterior median (95% CrI). 4,784 quarters, 1,196 cows.",
  "Covariates: front, left (within cow) and parity 2-7 (between cows); reference parity 1.",
  "",
  paste0("| | ", paste(names(fits), collapse = " | "), " |"),
  paste0("|---|", paste(rep("---", length(fits)), collapse = "|"), "|")
)
for (j in seq_along(cov_names_A)) {
  lines <- c(lines, paste0("| ", cov_names_A[j], " | ",
    paste(sapply(fits, par_row, par = "beta", idx = j), collapse = " | "), " |"))
}
lines <- c(lines,
  paste0("| theta | ", paste(sapply(fits, par_row, par = "theta"), collapse = " | "), " |"),
  paste0("| Var(log w) | ", paste(sapply(fits, par_row, par = "var_log_frailty"), collapse = " | "), " |"),
  paste0("| Kendall tau | ", paste(sapply(fits, par_row, par = "kendall_tau"), collapse = " | "), " |"),
  paste0("| Divergent transitions | ", paste(sapply(fits, ndiv), collapse = " | "), " |"),
  paste0("| Elapsed seconds | ", paste(sapply(fits, elap), collapse = " | "), " |"),
  "",
  paste0("cor(alpha, mean_log_w): ",
         paste(sprintf("%s %+.3f", names(rho), rho), collapse = " / ")),
  "",
  "## Within- vs between-cluster precision (spec B: front, left, multiparous)",
  "",
  "```", paste(capture.output(print(
    sB[, c("name", "scale", "median", "sd", "q5", "q95")])), collapse = "\n"), "```",
  "",
  sprintf("Mean posterior SD: within %.4f / between %.4f (ratio %.2f).", w_sd, b_sd, b_sd / w_sd),
  "For reference, category balance alone would give 1.06 and clustering alone 2.00."
)
writeLines(lines, "baselines.md")
cat("\n-> baselines.md written\n")
