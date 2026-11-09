# -----------------------------------------------------------------
# figures/fig1_baselines.R
# Figure 1: the three baseline constructions of Section 2.2.
#
# Illustrative curves, not fitted values. The figure makes one point:
# whether the cumulative baseline hazard is monotone by construction.
#   piecewise exponential : step hazard, piecewise-linear H0, never decreasing
#   I-spline              : non-negative coefficients on monotone basis functions
#   Royston-Parmar        : exp() guarantees positivity only; H0 can decrease
#
# For the Royston-Parmar panel, log H0 = a0 + b log t - a sin(w (log t - c)),
# and h0 < 0 somewhere requires b - a w < 0; with b = 1.45 and w = 2.1 this
# needs a > 0.690, and a = 0.90 is used. The three panels share the y-axis.
#
# Run from the figures/ directory: Rscript fig1_baselines.R
# Output: fig1_baselines.pdf (vector, 7 x 2.7 in)
# -----------------------------------------------------------------

OUT   <- "fig1_baselines.pdf"
tt    <- seq(0, 10, length.out = 400)
YMAX  <- 3

rescale <- function(y) y / max(y) * YMAX

# --- Piecewise exponential: cumulative of a step hazard ----------
pwe_cuts <- c(0, 2, 4, 6, 8, 10)
pwe_h    <- c(0.10, 0.22, 0.30, 0.52, 0.70)
pwe_H <- function(t) {
  H <- 0
  for (k in seq_along(pwe_h))
    if (t > pwe_cuts[k])
      H <- H + pwe_h[k] * (min(t, pwe_cuts[k + 1]) - pwe_cuts[k])
  H
}
H_pwe <- rescale(vapply(tt, pwe_H, 0))

# --- I-spline: smooth and monotone ------------------------------
H_isp <- rescale(0.035 * tt^1.6)

# --- Royston-Parmar: an example with a decreasing stretch -------
B_RP <- 1.45; W_RP <- 2.1; A_RP <- 0.90; C_RP <- log(4)
stopifnot(B_RP - A_RP * W_RP < 0)     # a decrease is possible
lt   <- log(pmax(tt, 1e-6))
H_rp <- exp(-2.3 + B_RP * lt - A_RP * sin(W_RP * (lt - C_RP)))
H_rp[1] <- 0
H_rp <- rescale(H_rp)

d      <- diff(H_rp)
dip    <- which(d < 0 & tt[-1] > 0.5)
dip_lo <- if (length(dip)) min(dip) else NA
dip_hi <- if (length(dip)) max(dip) + 1 else NA

# --- Drawing -----------------------------------------------------
pdf(OUT, width = 7, height = 2.7, pointsize = 9)
op <- par(mfrow = c(1, 3), mar = c(3.1, 3.3, 3.0, 0.8),
          mgp = c(2.0, 0.6, 0), las = 1, bty = "n")

panel <- function(y, main, sub, col) {
  plot(tt, y, type = "n", xlab = "t", ylab = expression(H[0](t)),
       ylim = c(0, YMAX * 1.08), xaxs = "i", yaxs = "i")
  lines(tt, y, col = col, lwd = 2.2)
  title(main = main, font.main = 1, line = 1.6)
  mtext(sub, side = 3, line = 0.4, cex = 0.72)
}

panel(H_pwe, "Piecewise exponential", "monotone by construction", "#0F6E56")
abline(v = pwe_cuts[2:5], col = "grey85", lty = 3)

panel(H_isp, "I-spline", "monotone by construction", "#534AB7")

panel(H_rp, "Royston-Parmar", "monotonicity not guaranteed", "#993C1D")
if (!is.na(dip_lo)) {
  lines(tt[dip_lo:dip_hi], H_rp[dip_lo:dip_hi], col = "#D85A30", lwd = 3.6)
  xm <- mean(tt[c(dip_lo, dip_hi)])
  ym <- min(H_rp[dip_lo:dip_hi])
  arrows(xm + 2.0, ym - 0.62, xm + 0.25, ym - 0.10,
         length = 0.05, col = "#D85A30")
  text(xm + 2.1, ym - 0.70, expression(h[0](t) < 0),
       col = "#D85A30", cex = 0.85, adj = 0)
}

par(op)
dev.off()

cat("->", OUT, "written\n")
if (!is.na(dip_lo)) {
  cat(sprintf("   decreasing stretch: t from %.2f to %.2f, H0 from %.2f to %.2f\n",
              tt[dip_lo], tt[dip_hi], H_rp[dip_lo], H_rp[dip_hi]))
} else {
  cat("   warning: no decreasing stretch; set A_RP above",
      sprintf("%.3f", B_RP / W_RP), "\n")
}
