# -----------------------------------------------------------------
# udder/01_data_weibull.R
# Udder quarter infection data (Goethals et al., 2009): data for the
# Weibull frailty models of Section 5 (Tables 2 and 3).
#
# Clusters are cows and the single covariate is multiparity, as in the
# published analysis, so that the replication in Table 2 uses the same
# data and the comparison of frailty distributions in Table 3 changes
# nothing but the frailty distribution.
#
# Time is converted from days to quarter-years (91.31 days), the scale of
# the published implementation. Censoring type (cens_type): 2 right,
# 4 interval. Quarters already infected at the first inspection
# (left == 0) are left-censored; since H0(0) = 0, the interval-censored
# contribution with L = 0 is exactly the left-censored one, so they need
# no separate type, as in the published implementation.
#
# Usage: source("01_data_weibull.R") from the udder/ directory.
# Creates d and stan_data.
# -----------------------------------------------------------------

library(UdderQuarterInfectionData)

data("udderquarterinfection")
d <- udderquarterinfection

# ---- Time scale of the published implementation: days -> quarter-years ----
SCALE <- 91.31            # = 365.25 / 4
d$L_use <- d$left  / SCALE
d$U_use <- d$right / SCALE

# ---- Censoring type ----
# status == 0: right-censored; status == 1: interval-censored (lower, upper]
d$cens_type <- ifelse(d$status == 1, 4L, 2L)

# For right-censored quarters the published implementation uses the upper
# column as the censoring time, so it is placed in L here.
d$L_stan <- ifelse(d$cens_type == 2, d$U_use, d$L_use)
d$U_stan <- d$U_use

# Interval-censored rows need L < U
stopifnot(all(d$L_stan[d$cens_type == 4] < d$U_stan[d$cens_type == 4]))

# ---- Cluster index (cluster = cow) ----
d$cluster_id <- as.integer(factor(d$cowid))
n_clusters   <- max(d$cluster_id)

# ---- Covariate: multiparity (lactation > 1) ----
d$x <- as.numeric(d$lactation > 1)

# ===================================================================
# Check the data structure
# ===================================================================
n_interval <- sum(d$cens_type == 4 & d$L_stan >  0)   # interval-censored
n_left     <- sum(d$cens_type == 4 & d$L_stan == 0)   # left-censored
n_right    <- sum(d$cens_type == 2)

cat("---- Data structure ----\n")
cat("quarters       :", nrow(d),     " (expected 4784)\n")
cat("cows           :", n_clusters,  " (expected 1196)\n")
cat("interval       :", n_interval,  " (expected 1811)\n")
cat("left-censored  :", n_left,      " (expected   68)\n")
cat("right-censored :", n_right,     " (expected 2905)\n")
cat("multiparous    :", round(mean(d$x), 4), "\n")

stopifnot(nrow(d) == 4784, n_clusters == 1196)
stopifnot(n_interval == 1811, n_left == 68, n_right == 2905)

cat("cluster sizes:\n"); print(table(table(d$cluster_id)))
cat("-> data structure as expected\n\n")

# ===================================================================
# Stan data, shared by the gamma and log-normal frailty models
# ===================================================================
stan_data <- list(
  N          = nrow(d),
  n_clusters = n_clusters,
  cluster    = d$cluster_id,
  cens_type  = d$cens_type,
  L          = d$L_stan,
  U          = d$U_stan,
  x          = d$x
)
