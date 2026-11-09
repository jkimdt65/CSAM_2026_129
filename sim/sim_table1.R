# -----------------------------------------------------------------
# sim/sim_table1.R   Simulation study of Section 4 (Table 1)
#
# The data of each replication are generated after set.seed(seed0 + rep),
# so results do not depend on how mclapply seeds its child processes and,
# within a configuration, the I-spline and piecewise exponential fits use
# the same datasets.
#
# Design (all configurations)
#   500 clusters, log-normal frailty with sigma_b = 0.8,
#   beta1 = beta2 = 0.5, switch time v ~ U(1, 4), 30% delayed entry,
#   Weibull baseline H0(t) = 0.02 t^1.5, about 32% right censoring,
#   150 replications, 4 chains x (700 warmup + 700 sampling),
#   adapt_delta = 0.9.
#   Fitted baselines: piecewise exponential (K = 6) and I-spline
#   (3 interior knots), as set in make_stan_data() of sim_functions.R.
#
# Configurations
#   s4k10 : cluster size 4, 10 inspections
#   s4k20 : cluster size 4, 20 inspections
#   s2k20 : cluster size 2, 20 inspections
#
# Each configuration x baseline is saved separately and skipped if its
# file already exists, so an interrupted run can be resumed.
# Trial run (3 replications, files prefixed test_):
#   SIM_TEST=1 Rscript sim_table1.R
# Full run:
#   setsid nohup Rscript sim_table1.R > sim_table1_log.txt 2>&1 < /dev/null &   (from sim/)
# Output: sim_table1_{s4k10,s4k20,s2k20}_{isp,pwe}.rds, sim_table1_results.md
# Software: CmdStan 2.37.0
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
library(splines2)
library(parallel)

# Settings and functions (gen_one, make_stan_data, ...)
source("sim_functions.R")

N_CLUSTER <- 500
N_REP_S   <- 150

TEST   <- Sys.getenv("SIM_TEST") == "1"
PREFIX <- if (TEST) "test_" else ""
if (TEST) N_REP_S <- 3

CFG <- list(
  s4k10 = list(label = "size 4, 10 inspections", size = 4, ninsp = 10, seed0 = 21000),
  s4k20 = list(label = "size 4, 20 inspections", size = 4, ninsp = 20, seed0 = 22000),
  s2k20 = list(label = "size 2, 20 inspections", size = 2, ninsp = 20, seed0 = 23000)
)

cat("Compiling Stan models\n")
mods <- list(pwe = cmdstan_model("../stan/sim_pwe_tdc.stan"),
             isp = cmdstan_model("../stan/sim_isp_tdc.stan"))

one <- function(rep_id, cf, baseline) {
  tmpd <- file.path(tempdir(), sprintf("t1_%s_%s_%d", cf, baseline, rep_id))
  dir.create(tmpd, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmpd, recursive = TRUE), add = TRUE)
  tryCatch({
    set.seed(CFG[[cf]]$seed0 + rep_id)
    d <- gen_one(N_CLUSTER, CFG[[cf]]$size, TRUE)
    sd_list <- make_stan_data(d, baseline, "full")
    fit <- mods[[baseline]]$sample(
      data = sd_list, chains = N_CHAINS, parallel_chains = 1,
      iter_warmup = 700, iter_sampling = 700,
      seed = CFG[[cf]]$seed0 + rep_id, adapt_delta = 0.9,
      refresh = 0, show_messages = FALSE, output_dir = tmpd)
    g <- function(p) {
      dm <- as.numeric(fit$draws(p, format = "draws_matrix")[, 1])
      c(median(dm), quantile(dm, .025), quantile(dm, .975))
    }
    b1 <- g("beta[1]"); b2 <- g("beta_tdc"); sg <- g("sigma_out")
    dg <- fit$diagnostic_summary()
    data.frame(cfg = cf, baseline = baseline, rep = rep_id,
               b1_med = b1[1], b1_lo = b1[2], b1_hi = b1[3],
               b2_med = b2[1], b2_lo = b2[2], b2_hi = b2[3],
               sg_med = sg[1], sg_lo = sg[2], sg_hi = sg[3],
               ndiv = sum(dg$num_divergent),
               maxtd = sum(dg$num_max_treedepth),
               censrate = mean(d$ct == 2), n_obs = nrow(d), row.names = NULL)
  }, error = function(e) {
    data.frame(cfg = cf, baseline = baseline, rep = rep_id, err = conditionMessage(e))
  })
}

# I-spline first (fast), then piecewise exponential
combos <- expand.grid(baseline = c("isp", "pwe"), cf = names(CFG),
                      stringsAsFactors = FALSE)
combos <- combos[order(combos$baseline != "isp"), ]

for (k in seq_len(nrow(combos))) {
  cf <- combos$cf[k]; bl <- combos$baseline[k]
  rds <- sprintf("%ssim_table1_%s_%s.rds", PREFIX, cf, bl)
  if (file.exists(rds)) { cat(format(Sys.time(), "%H:%M"), rds, "exists, skipped\n"); next }
  N_INSPECT <- CFG[[cf]]$ninsp        # read by gen_one() in the child processes
  cat(format(Sys.time(), "%H:%M"), "start:", CFG[[cf]]$label, "/", toupper(bl), "\n")
  t0 <- Sys.time()
  res <- mclapply(seq_len(N_REP_S), one, cf = cf, baseline = bl,
                  mc.cores = N_CORES, mc.preschedule = FALSE)
  errs <- unlist(lapply(res, function(r) if (is.data.frame(r) && !is.null(r$err)) r$err else
                                         if (!is.data.frame(r)) "no result from mclapply" else NULL))
  if (length(errs)) cat("  failed:", length(errs), " first message:", errs[1], "\n")
  res <- do.call(rbind, lapply(res, function(r) if (is.data.frame(r) && is.null(r$err)) r else NULL))
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  attr(res, "minutes") <- mins
  saveRDS(res, rds)
  cat(format(Sys.time(), "%H:%M"), "done:", nrow(res), "/", N_REP_S, "replications,",
      round(mins), "min\n")
}

# -----------------------------------------------------------------
# Summary in the layout of Table 1
# -----------------------------------------------------------------
summ <- function(med, lo, hi, truth) {
  c(bias = mean(med) - truth, cover = mean(lo <= truth & truth <= hi),
    rmse = sqrt(mean((med - truth)^2)), mcse = sd(med) / sqrt(length(med)))
}
cols <- list()
for (cf in names(CFG)) for (bl in c("isp", "pwe")) {
  R <- readRDS(sprintf("%ssim_table1_%s_%s.rds", PREFIX, cf, bl))
  cols[[paste(cf, bl)]] <- list(
    n = nrow(R), mins = attr(R, "minutes"),
    b1 = summ(R$b1_med, R$b1_lo, R$b1_hi, BETA1_TRUE),
    b2 = summ(R$b2_med, R$b2_lo, R$b2_hi, BETA2_TRUE),
    sg = summ(R$sg_med, R$sg_lo, R$sg_hi, SIGMA_TRUE),
    div = mean(R$ndiv > 0), td = mean(R$maxtd > 0),
    cens = mean(R$censrate), nobs = mean(R$n_obs))
}
row <- function(lab, f) paste0("| ", lab, " | ", paste(vapply(cols, f, ""), collapse = " | "), " |")
hdr <- toupper(names(cols))
L <- c("# Simulation results (Table 1)", "",
  "500 clusters, sigma_b = 0.8, beta1 = beta2 = 0.5, Weibull baseline (shape 1.5),",
  "150 replications. Bias and RMSE refer to the posterior median.", "",
  paste0("| | ", paste(hdr, collapse = " | "), " |"),
  paste0("|", paste(rep("---", length(cols) + 1), collapse = "|"), "|"),
  row("Replications", function(x) as.character(x$n)),
  row("Mean n obs", function(x) sprintf("%.0f", x$nobs)),
  row("Right censoring", function(x) sprintf("%.1f%%", 100 * x$cens)),
  row("Bias beta1", function(x) sprintf("%+.3f", x$b1["bias"])),
  row("Bias beta2", function(x) sprintf("%+.3f", x$b2["bias"])),
  row("Bias sigma_b", function(x) sprintf("%+.3f", x$sg["bias"])),
  row("Coverage beta1", function(x) sprintf("%.3f", x$b1["cover"])),
  row("Coverage beta2", function(x) sprintf("%.3f", x$b2["cover"])),
  row("Coverage sigma_b", function(x) sprintf("%.3f", x$sg["cover"])),
  row("RMSE beta1", function(x) sprintf("%.3f", x$b1["rmse"])),
  row("RMSE beta2", function(x) sprintf("%.3f", x$b2["rmse"])),
  row("RMSE sigma_b", function(x) sprintf("%.3f", x$sg["rmse"])),
  row("MCSE bias beta1", function(x) sprintf("%.3f", x$b1["mcse"])),
  row("MCSE bias beta2", function(x) sprintf("%.3f", x$b2["mcse"])),
  row("MCSE bias sigma_b", function(x) sprintf("%.3f", x$sg["mcse"])),
  row("Reps with any divergence", function(x) sprintf("%.1f%%", 100 * x$div)),
  row("Reps reaching max treedepth", function(x) sprintf("%.1f%%", 100 * x$td)),
  row("Minutes", function(x) sprintf("%.0f", x$mins)))
writeLines(L, sprintf("%ssim_table1_results.md", PREFIX))
cat(L, sep = "\n")
