# -----------------------------------------------------------------
# sim/sim_shapes.R   Additional simulation (Supplementary Table S4)
#
# The main configuration of Table 1 (500 clusters of 4, log-normal frailty
# with sigma_b = 0.8, beta1 = beta2 = 0.5, v ~ U(1, 4), 30% delayed entry,
# 10 inspections, 150 replications) with two other generating baselines,
# their scales chosen to keep right censoring near 32%:
#   dec : decreasing,  H0(t) = 0.1508 t^0.7               (Weibull, shape 0.7)
#   uni : unimodal,    H0(t) = 0.1733 log(1 + (t/3)^3)    (log-logistic, shape 3;
#                                                          hazard peaks near t = 3.8)
# Fitted baselines as in Table 1: piecewise exponential (K = 6) and I-spline
# (3 interior knots), from make_stan_data() of sim_functions.R.
#
# Data generation follows gen_one() of sim_functions.R with H0 and its
# inverse as arguments. Each replication is generated after
# set.seed(seed0 + rep). Each scenario x baseline is saved separately and
# skipped if its file exists, so an interrupted run can be resumed.
#
# Trial run (3 replications, files prefixed test_):
#   SIM_TEST=1 Rscript sim_shapes.R
# Full run (from sim/):
#   setsid nohup Rscript sim_shapes.R > sim_shapes_log.txt 2>&1 < /dev/null &
# Output: sim_shapes_{dec,uni}_{isp,pwe}.rds, sim_shapes_results.md
# Software: CmdStan 2.37.0
# -----------------------------------------------------------------

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
library(splines2)
library(parallel)

# Settings and functions (gen_one, make_stan_data, ...)
source("sim_functions.R")

N_CLUSTER <- 500
SZ        <- 4
N_INSPECT <- 10
N_REP_S   <- 150

# Trial run: SIM_TEST=1 -> 3 replications per combination, files prefixed test_
TEST   <- Sys.getenv("SIM_TEST") == "1"
PREFIX <- if (TEST) "test_" else ""
if (TEST) N_REP_S <- 3

SCEN <- list(
  dec = list(label = "decreasing (Weibull 0.7)", seed0 = 11000,
             H0 = function(t) 0.1508 * t^0.7,
             H0inv = function(H) (H / 0.1508)^(1 / 0.7)),
  uni = list(label = "unimodal (log-logistic 3)", seed0 = 12000,
             H0 = function(t) 0.1733 * log1p((t / 3)^3),
             H0inv = function(H) 3 * expm1(H / 0.1733)^(1 / 3))
)

# As gen_one() of sim_functions.R, with the baseline as an argument
gen_one_shape <- function(n_clust, size, H0, H0inv) {
  N <- n_clust * size
  cid <- rep(seq_len(n_clust), each = size)
  w <- rep(exp(SIGMA_TRUE * rnorm(n_clust)), each = size)
  x1 <- rbinom(N, 1, 0.5)
  vchg <- runif(N, 1, 4)

  U <- runif(N); target <- -log(U)
  base <- w * exp(BETA1_TRUE * x1)
  Hv <- H0(vchg) * base
  tt <- numeric(N)
  pre <- target < Hv
  tt[pre] <- H0inv(target[pre] / base[pre])
  post <- !pre
  rem <- (target[post] - Hv[post]) / (base[post] * exp(BETA2_TRUE))
  tt[post] <- H0inv(H0(vchg[post]) + rem)

  A <- rep(0, N)
  lt_idx <- sample(N, round(LT_FRAC * N))
  A[lt_idx] <- runif(length(lt_idx), 0.5, 2.5)
  keep <- tt > A

  insp <- sort(c(0, runif(N_INSPECT - 1, 0, MAX_T), MAX_T))
  L <- numeric(N); Uu <- numeric(N); ct <- integer(N)
  for (i in seq_len(N)) {
    grid <- insp[insp > A[i]]
    if (length(grid) == 0) { ct[i] <- 2L; L[i] <- A[i]; Uu[i] <- 1e4; next }
    if (tt[i] > max(grid)) { ct[i] <- 2L; L[i] <- max(grid); Uu[i] <- 1e4
    } else {
      k <- findInterval(tt[i], grid)
      lo <- if (k == 0) A[i] else grid[k]
      hi <- grid[k + 1]
      if (lo <= A[i] + 1e-9) {
        ct[i] <- if (A[i] > 0) 4L else 3L
        L[i]  <- A[i]
        Uu[i] <- hi
      } else { ct[i] <- 4L; L[i] <- lo; Uu[i] <- hi }
    }
  }
  data.frame(cid = cid[keep], x1 = x1[keep], vchg = vchg[keep],
             A = A[keep], L = L[keep], U = Uu[keep], ct = ct[keep])
}

cat("Compiling Stan models\n")
mods <- list(pwe = cmdstan_model("../stan/sim_pwe_tdc.stan"),
             isp = cmdstan_model("../stan/sim_isp_tdc.stan"))

one <- function(rep_id, sc, baseline) {
  tmpd <- file.path(tempdir(), sprintf("sh_%s_%s_%d", sc, baseline, rep_id))
  dir.create(tmpd, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmpd, recursive = TRUE), add = TRUE)
  tryCatch({
    set.seed(SCEN[[sc]]$seed0 + rep_id)
    d <- gen_one_shape(N_CLUSTER, SZ, SCEN[[sc]]$H0, SCEN[[sc]]$H0inv)
    sd_list <- make_stan_data(d, baseline, "full")
    fit <- mods[[baseline]]$sample(
      data = sd_list, chains = N_CHAINS, parallel_chains = 1,
      iter_warmup = 700, iter_sampling = 700,
      seed = SCEN[[sc]]$seed0 + rep_id, adapt_delta = 0.9,
      refresh = 0, show_messages = FALSE, output_dir = tmpd)
    g <- function(p) {
      dm <- as.numeric(fit$draws(p, format = "draws_matrix")[, 1])
      c(median(dm), quantile(dm, .025), quantile(dm, .975))
    }
    b1 <- g("beta[1]"); b2 <- g("beta_tdc"); sg <- g("sigma_out")
    data.frame(scen = sc, baseline = baseline, rep = rep_id,
               b1_med = b1[1], b1_lo = b1[2], b1_hi = b1[3],
               b2_med = b2[1], b2_lo = b2[2], b2_hi = b2[3],
               sg_med = sg[1], sg_lo = sg[2], sg_hi = sg[3],
               ndiv = sum(fit$diagnostic_summary()$num_divergent),
               maxtd = max(fit$diagnostic_summary()$num_max_treedepth),
               censrate = mean(d$ct == 2), n_obs = nrow(d), row.names = NULL)
  }, error = function(e) {
    data.frame(scen = sc, baseline = baseline, rep = rep_id, err = conditionMessage(e))
  })
}

# I-spline first (fast), then piecewise exponential
combos <- expand.grid(baseline = c("isp", "pwe"), sc = names(SCEN), stringsAsFactors = FALSE)
combos <- combos[order(combos$baseline != "isp"), ]
secs <- c()

for (k in seq_len(nrow(combos))) {
  sc <- combos$sc[k]; bl <- combos$baseline[k]
  rds <- sprintf("%ssim_shapes_%s_%s.rds", PREFIX, sc, bl)
  if (file.exists(rds)) { cat(format(Sys.time(), "%H:%M"), rds, "exists, skipped\n"); next }
  cat(format(Sys.time(), "%H:%M"), "start:", SCEN[[sc]]$label, "/", toupper(bl), "\n")
  t0 <- Sys.time()
  res <- mclapply(seq_len(N_REP_S), one, sc = sc, baseline = bl,
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

# ===================================================================
# Summary in the layout of Table 1
# ===================================================================
summ <- function(med, lo, hi, truth) {
  c(bias = mean(med) - truth, cover = mean(lo <= truth & truth <= hi),
    rmse = sqrt(mean((med - truth)^2)), mcse = sd(med) / sqrt(length(med)))
}
cols <- list()
for (sc in names(SCEN)) for (bl in c("isp", "pwe")) {
  R <- readRDS(sprintf("%ssim_shapes_%s_%s.rds", PREFIX, sc, bl))
  cols[[paste(sc, bl)]] <- list(
    n = nrow(R), mins = attr(R, "minutes"),
    b1 = summ(R$b1_med, R$b1_lo, R$b1_hi, BETA1_TRUE),
    b2 = summ(R$b2_med, R$b2_lo, R$b2_hi, BETA2_TRUE),
    sg = summ(R$sg_med, R$sg_lo, R$sg_hi, SIGMA_TRUE),
    div = mean(R$ndiv > 0), td = mean(R$maxtd > 0),
    cens = mean(R$censrate), nobs = mean(R$n_obs))
}
row <- function(lab, f) paste0("| ", lab, " | ", paste(vapply(cols, f, ""), collapse = " | "), " |")
L <- c("# Additional simulation: decreasing and unimodal generating baselines", "",
  "Design as Table 1 (500 clusters of 4, sigma_b = 0.8, 10 inspections, 150 replications);",
  "only the generating baseline differs. Bias and RMSE refer to the posterior median.", "",
  paste0("| | ", paste(c("Dec ISP", "Dec PWE", "Uni ISP", "Uni PWE"), collapse = " | "), " |"),
  "|---|---|---|---|---|",
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
  row("Reps hitting max treedepth", function(x) sprintf("%.1f%%", 100 * x$td)),
  row("Minutes", function(x) sprintf("%.0f", x$mins)))
writeLines(L, paste0(PREFIX, "sim_shapes_results.md"))
cat("\n->", paste0(PREFIX, "sim_shapes_results.md"), "written\n")
