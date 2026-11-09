# -----------------------------------------------------------------
# bcos/01_data.R
# Breast cosmesis data (Finkelstein and Wolfe, 1985): data preparation
# for the three baseline constructions of Section 2.2.
#
# The data are taken as bcos2 from the CRAN package survivalMPL. Its
# coding differs from interval::bcos: left-censored rows have left = NA
# and right-censored rows have right = NA. These are recoded below to
# L = 0 and U = BIG, as the Stan programs expect.
#
# Censoring type (cens_type): 1 exact, 2 right, 3 left, 4 interval.
# There are no exact observations in these data (checked below).
#
# Usage: source("01_data.R") from the bcos/ directory. Creates
#   d, X, BIG, cuts, K, ns_basis, ctr, pos_rp,
#   stan_data_pwe_bcos, stan_data_isp_bcos, stan_data_rp_bcos
# -----------------------------------------------------------------

library(survivalMPL)
library(splines2)
library(splines)

data(bcos2, package = "survivalMPL")
d <- bcos2

# ---- Check against the published description of the data ----
n_left  <- sum(is.na(d$left))
n_right <- sum(is.na(d$right))
n_int   <- sum(!is.na(d$left) & !is.na(d$right))
n_exact <- sum(!is.na(d$left) & !is.na(d$right) & d$left == d$right, na.rm = TRUE)

cat("---- Data structure ----\n")
cat("observations   :", nrow(d), "\n")
cat("left-censored  :", n_left, "(expected 5)\n")
cat("right-censored :", n_right, "(expected 38)\n")
cat("interval       :", n_int, "(expected 51)\n")
cat("exact          :", n_exact, "(expected 0)\n")
print(table(d$treatment))

stopifnot(nrow(d) == 94, n_left == 5, n_right == 38, n_int == 51, n_exact == 0)
stopifnot(sum(d$treatment == "Rad") == 46, sum(d$treatment == "RadChem") == 48)

# Time is in months, as in the original data; no rescaling is needed.

# ---- L, U and censoring type ----
BIG <- 1e4    # large finite value standing for +Inf (follow-up is 0 to 48 months)

d$cens_type <- ifelse(is.na(d$left), 3L,                # left-censored
               ifelse(is.na(d$right), 2L,                # right-censored
                                       4L))               # interval-censored

d$L_stan <- ifelse(d$cens_type == 3, 0,                  # left:  L = 0
            ifelse(d$cens_type == 2, d$left, d$left))    # right, interval: L = left
d$U_stan <- ifelse(d$cens_type == 2, BIG,                 # right: U = BIG
            ifelse(d$cens_type == 3, d$right, d$right))   # left, interval: U = right

stopifnot(all(!is.na(d$L_stan)), all(!is.na(d$U_stan)))
stopifnot(all(d$L_stan[d$cens_type == 4] < d$U_stan[d$cens_type == 4]))
stopifnot(all(d$L_stan >= 0), all(d$U_stan > d$L_stan | d$cens_type != 4))

cat("L range:", round(range(d$L_stan[d$cens_type != 2]), 2), "\n")
cat("U range (excluding right-censored):", round(range(d$U_stan[d$cens_type != 2]), 2), "\n")

# ---- Covariate: radiotherapy plus adjuvant chemotherapy ----
d$radchem <- as.integer(d$treatment == "RadChem")
X <- matrix(d$radchem, ncol = 1)
colnames(X) <- "radchem"

# ===================================================================
# Piecewise exponential cutpoints: equally spaced quantiles of the
# endpoints of the Turnbull intersection intervals (Section 2.2)
# ===================================================================
turnbull_intervals <- function(L, U) {
  Ls <- sort(unique(L))
  Us <- sort(unique(U[U < BIG]))
  cand <- expand.grid(l = Ls, u = Us)
  cand <- cand[cand$l < cand$u, ]
  keep <- vapply(seq_len(nrow(cand)), function(i) {
    l <- cand$l[i]; u <- cand$u[i]
    !any(Ls > l & Ls < u) && !any(Us > l & Us < u)
  }, logical(1))
  out <- cand[keep, ]; out <- out[order(out$l), ]; rownames(out) <- NULL
  out
}

make_pwe_cuts <- function(L, U, K) {
  tb <- turnbull_intervals(L, U)
  ep <- sort(unique(c(tb$l, tb$u))); ep <- ep[ep > 0 & ep < BIG]
  inner <- quantile(ep, probs = seq(1, K - 1) / K)
  upper <- max(ep) * 1.05
  cuts <- sort(unique(c(0, unname(inner), upper)))
  cuts
}

K <- 6   # fewer intervals than for the udder data (K = 8): 94 patients against 4,784 quarters
cuts <- make_pwe_cuts(d$L_stan, d$U_stan, K)
K <- length(cuts) - 1

cat("\n---- Piecewise exponential cutpoints (K =", K, ") ----\n")
print(round(cuts, 3))
ev <- d$U_stan[d$cens_type == 4]
cat("events per interval (by upper endpoint):\n")
print(table(cut(ev, breaks = cuts, include.lowest = TRUE)))

# ===================================================================
# I-spline basis: knots at quantiles of the observed interval endpoints
# ===================================================================
pos <- c(d$L_stan[d$L_stan > 0], d$U_stan[d$U_stan < BIG])
bnd <- c(0, max(pos) * 1.05)
nk  <- 3   # fewer interior knots than for the udder data (5)
ik  <- quantile(pos, probs = seq_len(nk) / (nk + 1))
ib  <- iSpline(d$L_stan, knots = ik, degree = 3, intercept = TRUE,
               Boundary.knots = bnd)
mb  <- mSpline(d$L_stan, knots = ik, degree = 3, intercept = TRUE,
               Boundary.knots = bnd)   # M-spline for h0(t); used only for exact observations
IL <- predict(ib, d$L_stan)
IU <- predict(ib, pmin(d$U_stan, bnd[2]))   # right-censored rows: U = BIG -> boundary
ML <- predict(mb, d$L_stan)                 # not referenced here (no exact observations)
stopifnot(all(abs(IL[d$L_stan == 0, ]) < 1e-10))
nbasis <- ncol(IL)
cat("\nI-spline basis functions:", nbasis, "\n")

# ===================================================================
# Royston-Parmar basis: natural cubic spline in log t, centred so that
# the intercept alpha carries the level of log H0 (Section 2.2)
# ===================================================================
pos_rp  <- c(d$L_stan[d$L_stan > 0], d$U_stan[d$U_stan < BIG])
log_pos <- log(pos_rp)
df_rp   <- 3   # fewer degrees of freedom than for the udder data (5)
ns_basis <- ns(log_pos, df = df_rp, intercept = FALSE)
ctr <- colMeans(predict(ns_basis, log_pos))

L_zero <- as.integer(d$L_stan == 0)
NS_L <- matrix(0, nrow(d), ncol(ns_basis))
NS_L[L_zero == 0, ] <- sweep(predict(ns_basis, log(d$L_stan[L_zero == 0])), 2, ctr, "-")
NS_U <- sweep(predict(ns_basis, log(pmin(d$U_stan, max(pos_rp)))), 2, ctr, "-")
NS_ent <- matrix(0, nrow(d), ncol(ns_basis))   # no delayed entry in these data
dNS_L  <- matrix(0, nrow(d), ncol(ns_basis))   # derivative basis; used only for exact observations

# Rows with L = 0 are all left-censored, and the left-censored branch of the
# Stan program uses only H0(U), so NS_L is never evaluated at log(0).
stopifnot(all(d$L_stan[d$cens_type == 4] > 0))
stopifnot(all(NS_L[L_zero == 1, ] == 0))

cat("Royston-Parmar degrees of freedom:", ncol(ns_basis), "\n")

# ===================================================================
# Stan data
# ===================================================================
base_list <- list(N = nrow(d), P = ncol(X), X = X,
                  cens_type = d$cens_type, L = d$L_stan, U = d$U_stan,
                  entry = rep(0, nrow(d)))

stan_data_pwe_bcos <- c(base_list, list(K = K, cuts = cuts))
stan_data_isp_bcos <- c(base_list, list(nbasis = nbasis, IL = IL, IU = IU,
                                        Ient = matrix(0, nrow(d), nbasis),
                                        ML = ML))
stan_data_rp_bcos  <- c(base_list, list(df = df_rp, NS_L = NS_L, NS_U = NS_U,
                                        NS_ent = NS_ent, dNS_L = dNS_L))

cat("\n-> stan_data_pwe_bcos, stan_data_isp_bcos, stan_data_rp_bcos ready\n")
