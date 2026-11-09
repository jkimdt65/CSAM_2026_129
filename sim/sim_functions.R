# -----------------------------------------------------------------
# sim/sim_functions.R
# Settings and functions shared by the simulation scripts of Section 4
# (sim_table1.R, Table 1; sim_shapes.R, Supplementary Table S4).
#
#   gen_one()        : one simulated dataset (Weibull baseline, exogenous
#                      time-varying covariate, delayed entry, log-normal
#                      frailty, interval censoring on a shared grid)
#   make_stan_data() : Stan data for the piecewise exponential (K = 6)
#                      or I-spline (3 interior knots) fit
#
# True values: beta1 = 0.5 (time-fixed), beta2 = 0.5 (time-varying),
# sigma_b = 0.8. N_CLUSTER and N_INSPECT are reset by the calling script.
#
# Usage: source("sim_functions.R") from the sim/ directory.
# -----------------------------------------------------------------

library(cmdstanr)
library(splines2)
library(parallel)

set.seed(2026)

# ===================================================================
# Settings
# ===================================================================
N_CLUSTER    <- 200     # reset by the calling script
SIGMA_TRUE   <- 0.8     # log-normal frailty SD; Var(log w) = 0.64
BETA1_TRUE   <- 0.5     # time-fixed covariate
BETA2_TRUE   <- 0.5     # time-varying covariate
WEIB_SHAPE   <- 1.5
WEIB_SCALE   <- 0.02    # H0(t) = 0.02 t^1.5; about 32% right censoring
LT_FRAC      <- 0.30    # share of subjects with delayed entry
N_INSPECT    <- 6       # number of inspections; reset by the calling script
MAX_T        <- 10
N_CORES      <- 12
N_CHAINS     <- 4

# ===================================================================
# Data generation
#   time-varying covariate: x2 switches from 0 to 1 at time v (exogenous)
#   delayed entry: some subjects enter at A > 0; events before A are never
#   observed, so those subjects are not in the sample
# ===================================================================
gen_one <- function(n_clust, size, with_frailty) {
  N <- n_clust * size
  cid <- rep(seq_len(n_clust), each = size)

  # Log-normal frailty, log w = sigma_b * z with z ~ N(0, 1): the same form
  # as the non-centred parameterisation of the fitted model.
  w <- if (with_frailty) {
    rep(exp(SIGMA_TRUE * rnorm(n_clust)), each = size)
  } else rep(1, N)

  x1 <- rbinom(N, 1, 0.5)                    # time-fixed
  vchg <- runif(N, 1, 4)                     # time at which x2 switches to 1

  H0 <- function(t) WEIB_SCALE * t^WEIB_SHAPE
  H0inv <- function(H) (H / WEIB_SCALE)^(1/WEIB_SHAPE)

  U <- runif(N); target <- -log(U)
  base <- w * exp(BETA1_TRUE * x1)
  Hv <- H0(vchg) * base                       # cumulative hazard up to v (x2 = 0)
  tt <- numeric(N)
  pre <- target < Hv
  tt[pre] <- H0inv(target[pre] / base[pre])
  post <- !pre
  rem <- (target[post] - Hv[post]) / (base[post] * exp(BETA2_TRUE))
  tt[post] <- H0inv(H0(vchg[post]) + rem)

  # Delayed entry: A > 0 for a share of subjects; events before A are dropped
  A <- rep(0, N)
  lt_idx <- sample(N, round(LT_FRAC * N))
  A[lt_idx] <- runif(length(lt_idx), 0.5, 2.5)
  keep <- tt > A

  # Interval censoring on an inspection grid shared within the dataset
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
        # Event before the first inspection after entry: left-censored if
        # A = 0, interval-censored on (A, hi] if A > 0. Setting L = 0 in the
        # latter case would discard the condition T > A, so L = A.
        ct[i] <- if (A[i] > 0) 4L else 3L
        L[i]  <- A[i]
        Uu[i] <- hi
      } else { ct[i] <- 4L; L[i] <- lo; Uu[i] <- hi }
    }
  }

  data.frame(cid = cid[keep], x1 = x1[keep], vchg = vchg[keep],
             A = A[keep], L = L[keep], U = Uu[keep], ct = ct[keep])
}

# ===================================================================
# Stan data
# ===================================================================
make_stan_data <- function(d, baseline, tdc = "full") {
  N <- nrow(d)
  cid <- as.integer(factor(d$cid))
  X <- cbind(x1 = d$x1)

  # tdc = "base" moves the switch beyond the end of follow-up so that x2
  # stays 0 (used only in exploratory runs). The piecewise exponential
  # model uses vchg directly and the I-spline model uses Ivchg and the
  # after_v flags, so both read vv.
  vv <- if (tdc == "base") rep(1e4, N) else d$vchg

  # The fitted model always includes the frailty.
  base <- list(N = N, P = ncol(X), X = X, cens_type = d$ct,
               L = d$L, U = d$U, entry = d$A,
               n_clusters = max(cid), cluster = cid,
               vchg = vv, use_frailty = 1L)

  if (baseline == "pwe") {
    ep <- sort(unique(c(d$L[d$L > 0], d$U[d$U < 1e4])))
    K <- 6
    cuts <- sort(unique(c(0, quantile(ep, seq(1, K-1)/K), max(ep) * 1.05)))
    base <- c(base, list(K = length(cuts) - 1, cuts = cuts))
  } else {
    # Knots at quartiles of the positive observed endpoints and the switch
    # times; boundary knots from 0 to just beyond the largest of these.
    pos <- c(d$L[d$L > 0], d$U[d$U < 1e4], d$vchg)
    bnd <- c(0, max(pos) * 1.05)
    ik  <- quantile(pos, probs = (1:3)/4)
    ib  <- iSpline(d$L, knots = ik, degree = 3, intercept = TRUE,
                   Boundary.knots = bnd)
    Ucap <- pmin(d$U, bnd[2])
    vcap <- pmin(vv, bnd[2])
    base <- c(base, list(
      nbasis = ncol(predict(ib, d$L)),
      IL    = predict(ib, d$L),
      IU    = predict(ib, Ucap),
      Ivchg = predict(ib, vcap),
      Ient  = predict(ib, pmin(d$A, bnd[2])),
      L_after_v   = as.integer(d$L  > vv),
      U_after_v   = as.integer(Ucap > vv),
      ent_after_v = as.integer(d$A  > vv)))
  }
  base
}
