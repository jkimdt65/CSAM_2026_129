# -----------------------------------------------------------------
# udder/05_tables_95.R
# Tables 2, 3 and 4 of the paper, with 95% credible intervals, from the
# fits saved by 02_fit_frailty.R and 04_fit_baselines.R. No refitting.
#
# 95% intervals are used throughout so that the posterior intervals of
# Table 2 can be set beside the 95% confidence intervals of the published
# maximum likelihood fit.
#
# Run from the udder/ directory, after 02 and 04:
#   Rscript 05_tables_95.R
# Output: tables_95.md
# -----------------------------------------------------------------

library(cmdstanr)
library(loo)

P <- 0.95
lo <- (1 - P) / 2; hi <- 1 - lo

need <- c("fit_gamma.rds", "fit_lognormal.rds", "loo_frailty.rds",
          "fit_wb_A.rds", "fit_pwe_A.rds", "fit_isp_A.rds", "fit_wb_B.rds")
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing files: ", paste(miss, collapse = ", "))

get <- function(fit, par, idx = 1) {
  as.numeric(fit$draws(par, format = "draws_matrix")[, idx])
}
f3 <- function(v) sprintf("%.3f (%.3f, %.3f)",
                          median(v), quantile(v, lo), quantile(v, hi))
ndiv <- function(fit) sum(fit$diagnostic_summary()$num_divergent)

FG <- readRDS("fit_gamma.rds"); FL <- readRDS("fit_lognormal.rds")
LO <- readRDS("loo_frailty.rds")
WA <- readRDS("fit_wb_A.rds"); PA <- readRDS("fit_pwe_A.rds")
IA <- readRDS("fit_isp_A.rds"); WB <- readRDS("fit_wb_B.rds")

cov_A <- c("front","left","parity2","parity3","parity4","parity5","parity6","parity7")
cov_B <- c("front","left","multiparous")

# Published maximum likelihood fit (Goethals et al., 2009): estimates and
# variances, with Wald 95% intervals. The reference implementation took
# 1,038 seconds on one core for this fit.
mle <- data.frame(
  par = c("theta", "lambda", "gamma_shape", "beta"),
  est = c(3.7967246, 0.1201593, 1.9672298, 0.8590531),
  var = c(0.049281911, 0.0001982213, 0.0017306214, 0.020269244))
mle$lo <- mle$est - 1.96 * sqrt(mle$var)
mle$hi <- mle$est + 1.96 * sqrt(mle$var)
mrow <- function(p) with(mle[mle$par == p, ], sprintf("%.3f | (%.3f, %.3f)", est, lo, hi))

L <- c(sprintf("# Tables at %.0f%% credible intervals", 100 * P), "",
       "Posterior median (95% CrI), from the saved fits.", "")

# ---- Table 2 ----
L <- c(L, "## Table 2. Analytic marginalisation (published MLE) vs sampled frailties", "",
  "| Parameter | MLE | 95% CI | Posterior median (95% CrI) |",
  "|---|---|---|---|",
  paste0("| theta | ",        mrow("theta"),       " | ", f3(get(FG,"theta")), " |"),
  paste0("| lambda | ",       mrow("lambda"),      " | ", f3(get(FG,"lambda")), " |"),
  paste0("| Weibull shape | ", mrow("gamma_shape"), " | ", f3(get(FG,"gamma_shape")), " |"),
  paste0("| beta (multiparous) | ", mrow("beta"), " | ", f3(get(FG,"beta")), " |"),
  paste0("| Kendall tau | ", sprintf("%.3f", 3.7967246 / (3.7967246 + 2)), " | | ",
         f3(get(FG,"kendall_tau")), " |"),
  "", sprintf("Divergent transitions: %d. Stan elapsed: %.0f seconds on 4 cores.",
              ndiv(FG), LO$seconds[["gamma"]]), "")

# ---- Table 3 ----
el <- function(x) sprintf("%.1f (%.1f)", x$estimates["elpd_loo", "Estimate"],
                          x$estimates["elpd_loo", "SE"])
dl <- as.data.frame(LO$cmp_obs)
# loo >= 2.8 returns the model names in a column, older versions as row names
if ("model" %in% names(dl)) rownames(dl) <- dl$model
L <- c(L, "## Table 3. Gamma vs log-normal frailty (Weibull baseline)", "",
  "| | Gamma | Log-normal |", "|---|---|---|",
  paste0("| beta (multiparous) | ", f3(get(FG,"beta")), " | ", f3(get(FL,"beta")), " |"),
  paste0("| Hazard ratio | ", f3(get(FG,"hazard_ratio")), " | ", f3(get(FL,"hazard_ratio")), " |"),
  paste0("| Frailty parameter | theta = ", f3(get(FG,"theta")),
         " | sigma_b = ", f3(get(FL,"sigma_b")), " |"),
  paste0("| Var(log w) | ", f3(trigamma(1/get(FG,"theta"))),
         " | ", f3(get(FL,"var_log_frailty")), " |"),
  paste0("| Kendall tau | ", f3(get(FG,"kendall_tau")), " | ", f3(LO$tau_l), " |"),
  paste0("| elpd_loo (SE) | ", el(LO$loo_g), " | ", el(LO$loo_l), " |"),
  paste0("| elpd difference (SE) | | ",
         sprintf("%.1f (%.1f)", dl["lognormal", "elpd_diff"], dl["lognormal", "se_diff"]), " |"),
  paste0("| Pareto k > 0.7 | ", sum(pareto_k_values(LO$loo_g) > 0.7), " | ",
         sum(pareto_k_values(LO$loo_l) > 0.7), " |"),
  paste0("| Divergent transitions | ", ndiv(FG), " | ", ndiv(FL), " |"), "")

# ---- Table 4 ----
hdr <- c("| | Weibull | Piecewise exponential (K = 8) | I-spline |", "|---|---|---|---|")
rows <- character(0)
for (j in seq_along(cov_A))
  rows <- c(rows, sprintf("| %s | %s | %s | %s |", cov_A[j],
    f3(get(WA,"beta",j)), f3(get(PA,"beta",j)), f3(get(IA,"beta",j))))
rows <- c(rows,
  sprintf("| theta | %s | %s | %s |", f3(get(WA,"theta")), f3(get(PA,"theta")), f3(get(IA,"theta"))),
  sprintf("| Var(log w) | %s | %s | %s |",
          f3(get(WA,"var_log_frailty")), f3(get(PA,"var_log_frailty")), f3(get(IA,"var_log_frailty"))),
  sprintf("| Kendall tau | %s | %s | %s |",
          f3(get(WA,"kendall_tau")), f3(get(PA,"kendall_tau")), f3(get(IA,"kendall_tau"))),
  sprintf("| Divergent transitions | %d | %d | %d |", ndiv(WA), ndiv(PA), ndiv(IA)))
L <- c(L, "## Table 4. Baseline construction (gamma frailty, spec A)", "", hdr, rows, "")

# ---- spec B ----
sB <- WB$summary("beta"); sB$name <- cov_B
L <- c(L, "## Within- vs between-cluster precision (spec B, Weibull)", "",
  "| | Scale | Posterior median | 95% CrI | Posterior SD | Bulk ESS |",
  "|---|---|---|---|---|---|")
for (j in seq_along(cov_B)) {
  v <- get(WB, "beta", j)
  L <- c(L, sprintf("| %s | %s | %.3f | (%.3f, %.3f) | %.4f | %.0f |",
    cov_B[j], if (cov_B[j] %in% c("front","left")) "within" else "between",
    median(v), quantile(v, lo), quantile(v, hi), sd(v), sB$ess_bulk[j]))
}

# ---- Convergence ----
mrh <- function(fit) max(fit$summary()$rhat, na.rm = TRUE)
L <- c(L, "", "## Largest R-hat", "",
  sprintf("gamma %.4f, log-normal %.4f, Weibull A %.4f, piecewise exponential A %.4f, I-spline A %.4f, Weibull B %.4f",
          mrh(FG), mrh(FL), mrh(WA), mrh(PA), mrh(IA), mrh(WB)))

writeLines(L, "tables_95.md")
cat("-> tables_95.md written\n")
