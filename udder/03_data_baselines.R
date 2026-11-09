# -----------------------------------------------------------------
# udder/03_data_baselines.R
# Udder quarter infection data: data for the comparison of baseline
# constructions (Table 4, Supplementary Table S1).
#
# Time and censoring are coded exactly as in 01_data_weibull.R; only the
# covariates differ.
#   spec A: front, left (within cow) + parity dummies 2-7 (between cows)
#           -> Table 4
#   spec B: front, left (within cow) + multiparous (between cows)
#           -> within- vs between-cluster precision (Section 5)
#
# Piecewise exponential cutpoints are equally spaced quantiles of the
# endpoints of the Turnbull intersection intervals (K = 8); I-spline knots
# are quantiles of the positive observed endpoints (5 interior knots).
#
# Usage: source("03_data_baselines.R") from the udder/ directory. Creates
#   d, cov_names_A, cov_names_B, cuts, K, turnbull_intervals(),
#   make_pwe_cuts(), stan_data_wb_A, stan_data_wb_B, stan_data_pwe_A,
#   stan_data_isp_A
# -----------------------------------------------------------------

library(UdderQuarterInfectionData)
library(splines2)

# ===================================================================
# 0. Turnbull intersection intervals and piecewise exponential cutpoints
# ===================================================================
turnbull_intervals <- function(L, U) {
  Ls <- sort(unique(L))
  Us <- sort(unique(U[is.finite(U)]))

  cand <- expand.grid(l = Ls, u = Us)
  cand <- cand[cand$l < cand$u, ]

  keep <- vapply(seq_len(nrow(cand)), function(i) {
    l <- cand$l[i]; u <- cand$u[i]
    !any(Ls > l & Ls < u) && !any(Us > l & Us < u)
  }, logical(1))

  out <- cand[keep, ]
  out <- out[order(out$l), ]
  rownames(out) <- NULL
  out
}

make_pwe_cuts <- function(L, U, K, method = c("turnbull", "quantile")) {
  method <- match.arg(method)

  if (method == "quantile") {
    finite_times <- c(L, U[is.finite(U)])
    finite_times <- finite_times[finite_times > 0]
    inner <- quantile(finite_times, probs = seq(1, K - 1) / K)
    upper <- max(finite_times) * 1.05
  } else {
    tb <- turnbull_intervals(L, U)
    endpoints <- sort(unique(c(tb$l, tb$u)))
    endpoints <- endpoints[endpoints > 0]
    if (length(endpoints) < K) {
      warning("fewer Turnbull endpoints than K; using quantiles of the observed times")
      return(make_pwe_cuts(L, U, K, method = "quantile"))
    }
    inner <- quantile(endpoints, probs = seq(1, K - 1) / K)
    upper <- max(endpoints) * 1.05
  }

  cuts <- sort(unique(c(0, unname(inner), upper)))
  if (length(cuts) < K + 1) {
    warning("duplicate cutpoints; number of intervals reduced to ", length(cuts) - 1)
  }
  cuts
}

# ===================================================================
# 1. Time and censoring, as in 01_data_weibull.R
# ===================================================================
data("udderquarterinfection")
d <- udderquarterinfection

SCALE <- 91.31
d$L_use <- d$left  / SCALE
d$U_use <- d$right / SCALE

d$cens_type <- ifelse(d$status == 1, 4L, 2L)
d$L_stan    <- ifelse(d$cens_type == 2, d$U_use, d$L_use)
d$U_stan    <- d$U_use

stopifnot(all(d$L_stan[d$cens_type == 4] < d$U_stan[d$cens_type == 4]))

d$cluster_id <- as.integer(factor(d$cowid))
n_clusters   <- max(d$cluster_id)

cat("---- Data structure ----\n")
cat("quarters", nrow(d), "/ cows", n_clusters, "\n")
stopifnot(nrow(d) == 4784, n_clusters == 1196)

# ===================================================================
# 2. Covariates
# ===================================================================
# --- within cow (vary between quarters) ---
# quarter codes: LF left front, LR left rear, RF right front, RR right rear
d$front <- as.integer(d$quarter %in% c("LF", "RF"))
d$left  <- as.integer(d$quarter %in% c("LF", "LR"))

# --- between cows (constant within cow) ---
PARITY_LEVELS <- 2:7
for (p in PARITY_LEVELS) {
  d[[paste0("parity", p)]] <- as.integer(d$lactation == p)
}
d$multiparous <- as.integer(d$lactation > 1)   # spec B

cov_names_A <- c("front", "left", paste0("parity", PARITY_LEVELS))
cov_names_B <- c("front", "left", "multiparous")

X_A <- as.matrix(d[, cov_names_A]); storage.mode(X_A) <- "double"
X_B <- as.matrix(d[, cov_names_B]); storage.mode(X_B) <- "double"

stopifnot(ncol(X_A) == 8, ncol(X_B) == 3, !any(is.na(X_A)))
stopifnot(qr(X_A)$rank == 8, qr(X_B)$rank == 3)

cat("\n---- Covariates ----\n")
cat("front x left (should be balanced):\n")
print(table(front = d$front, left = d$left))

cat("\nquarters and cows by parity:\n")
tab <- table(d$lactation)
print(rbind(quarters = as.integer(tab), cows = as.integer(tab) / 4))
stopifnot(all(tab %% 4 == 0))

# front varies within every cow; parity is constant within every cow
within_var  <- tapply(d$front,   d$cluster_id, function(v) length(unique(v)))
between_var <- tapply(d$parity2, d$cluster_id, function(v) length(unique(v)))
cat("\nshare of cows in which front varies   :", mean(within_var > 1),  "(should be 1)\n")
cat("share of cows in which parity2 is constant:", mean(between_var == 1), "(should be 1)\n")
stopifnot(all(within_var > 1), all(between_var == 1))

# Balance of the binary covariates. Posterior SDs depend on how rare a
# category is as well as on whether it varies within or between cows; the
# parity dummies are very unbalanced (parity 7: 88 quarters), which is why
# the within/between comparison uses spec B.
cat("\n---- Covariate balance ----\n")
bal <- sapply(c(cov_names_A, "multiparous"), function(v) sum(d[[v]]))
print(rbind(n_ones = bal, prop = round(bal / nrow(d), 3)))

# ===================================================================
# 3. Piecewise exponential cutpoints
# ===================================================================
K <- 8
cuts <- make_pwe_cuts(L = d$L_stan, U = d$U_stan, K = K, method = "turnbull")
K    <- length(cuts) - 1

cat("\n---- Piecewise exponential cutpoints (K =", K, ") ----\n")
print(round(cuts, 4))

ev <- d$U_stan[d$cens_type == 4]
cat("events per interval (by upper endpoint):\n")
print(table(cut(ev, breaks = cuts, include.lowest = TRUE)))
# An interval without events would be informed by the random-walk prior
# alone; K should then be reduced.

# ===================================================================
# 4. I-spline basis
# ===================================================================
# Knots at equally spaced quantiles of all positive observed times (L > 0
# and U); the same basis is evaluated at L and at U, so that H0(L) and
# H0(U) are the same function. Boundary: 0 on the left, so that H0(0) = 0
# and rows with L = 0 need no special treatment, and 1.05 times the
# largest time on the right.
DEGREE <- 3      # cubic
NKNOTS <- 5      # interior knots

pos_times <- c(d$L_stan[d$L_stan > 0], d$U_stan)
boundary  <- c(0, max(pos_times) * 1.05)
iknots    <- quantile(pos_times, probs = seq_len(NKNOTS) / (NKNOTS + 1))

ibasis <- iSpline(d$L_stan, knots = iknots, degree = DEGREE,
                  intercept = TRUE, Boundary.knots = boundary)

IL   <- predict(ibasis, d$L_stan)
IU   <- predict(ibasis, d$U_stan)
Ient <- matrix(0, nrow(d), ncol(IL))     # no delayed entry: I(0) = 0
nbasis <- ncol(IL)

cat("\n---- I-spline basis ----\n")
cat("degree", DEGREE, "/ interior knots", NKNOTS, "-> nbasis", nbasis, "\n")
cat("interior knots (quarter-years):", round(iknots, 4), "\n")

stopifnot(all(abs(IL[d$L_stan == 0, ]) < 1e-12))           # H0(0) = 0
stopifnot(all(IL >= -1e-12), all(IL <= 1 + 1e-12))
stopifnot(all(IU >= -1e-12), all(IU <= 1 + 1e-12))
ic <- d$cens_type == 4
stopifnot(all(IU[ic, ] - IL[ic, ] >= -1e-10))                # I(U) >= I(L)
stopifnot(!any(is.na(IL)), !any(is.na(IU)))

# ===================================================================
# 5. Stan data
# ===================================================================
base_list <- list(
  N = nrow(d), n_clusters = n_clusters, cluster = d$cluster_id,
  cens_type = d$cens_type, L = d$L_stan, U = d$U_stan
)

stan_data_wb_A  <- c(base_list, list(P = ncol(X_A), X = X_A))
stan_data_wb_B  <- c(base_list, list(P = ncol(X_B), X = X_B))
stan_data_pwe_A <- c(stan_data_wb_A,
                     list(entry = rep(0, nrow(d)), K = K, cuts = cuts))
stan_data_isp_A <- c(stan_data_wb_A,
                     list(entry = rep(0, nrow(d)), nbasis = nbasis,
                          IL = IL, IU = IU, Ient = Ient))

cat("\n-> spec A (P = 8), spec B (P = 3), piecewise exponential (K = ", K,
    "), I-spline (nbasis = ", nbasis, ") ready\n", sep = "")
