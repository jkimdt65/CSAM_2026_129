# -----------------------------------------------------------------
# udder/fig2_frailty.R
# Figure 2: posterior medians and 95% credible intervals of the cow-level
# frailties from the gamma frailty fit (fit_gamma.rds), ordered by
# posterior median and coloured by the number of infected quarters.
#
# Vertical axis: for cows without an infected quarter the lower limits
# reach about 1e-9; plotting them all would compress the range 1e-2 to
# 1e+1 that carries the point, so the axis is cut at Y_FLOOR = 1e-4. The
# number of cut intervals and the smallest lower limit are printed for
# the caption.
#
# Order check: the frailties w[i] must be in the order of
# as.integer(factor(cowid)); if they are, the posterior medians correlate
# strongly with the number of infected quarters. The script stops if the
# rank correlation is below MIN_COR.
#
# Run from the udder/ directory, after 02_fit_frailty.R:
#   Rscript fig2_frailty.R
# Output: fig2_frailty.pdf (vector, 7 x 3.4 in)
# -----------------------------------------------------------------

library(cmdstanr)
library(UdderQuarterInfectionData)

FIT  <- "fit_gamma.rds"
OUT  <- "fig2_frailty.pdf"
PROB <- c(0.025, 0.975)          # 95% credible intervals
MIN_COR <- 0.30                  # threshold for the order check
Y_FLOOR <- 1e-4                  # lower limit of the vertical axis

# --- Cow index, built as in 01_data_weibull.R ---
data("udderquarterinfection")
d <- udderquarterinfection
d$cluster_id <- as.integer(factor(d$cowid))
d$infected   <- as.integer(d$status == 1)

stopifnot(nrow(d) == 4784, max(d$cluster_id) == 1196)

n_inf <- as.integer(tapply(d$infected, d$cluster_id, sum))   # 0..4 per cow

# --- Frailties from the fit ---
fit <- readRDS(FIT)
W   <- fit$draws("w", format = "draws_matrix")
stopifnot(ncol(W) == length(n_inf))

med <- apply(W, 2, median)
lo  <- apply(W, 2, quantile, probs = PROB[1])
hi  <- apply(W, 2, quantile, probs = PROB[2])

# --- Order check ---
rho <- cor(med, n_inf, method = "spearman")
cat(sprintf("rank correlation of posterior median and infected quarters: %.3f\n", rho))
if (is.na(rho) || rho < MIN_COR) {
  stop(sprintf(paste0("rank correlation %.3f is low: the order of w in ",
                      "fit_gamma.rds may differ from as.integer(factor(cowid))."), rho))
}

# --- Ordering and colours ---
ord <- order(med)
med <- med[ord]; lo <- lo[ord]; hi <- hi[ord]; k <- n_inf[ord]

pal <- c("#C9D4D0", "#8FB8AC", "#4E9080", "#22685A", "#0B4034")  # 0..4
col <- pal[k + 1]

# --- Drawing ---
pdf(OUT, width = 7, height = 3.4, pointsize = 9)
op <- par(mar = c(3.4, 3.8, 0.8, 0.8), mgp = c(2.4, 0.6, 0), las = 1,
          bty = "n")

x  <- seq_along(med)
yl <- range(c(pmax(lo, Y_FLOOR), hi))

plot(x, med, type = "n", log = "y", ylim = yl,
     xlab = "Cow, ordered by posterior median frailty",
     ylab = expression(w[i]))
segments(x, pmax(lo, Y_FLOOR), x, hi, col = col, lwd = 0.5)
lines(x, med, col = "grey20", lwd = 0.9)
abline(h = 1, col = "grey45", lty = 2)

legend("topleft", bty = "n", cex = 0.82, horiz = TRUE,
       legend = 0:4, fill = pal, border = NA,
       title = "Infected quarters", title.adj = 0)

par(op)
dev.off()

cat("->", OUT, "written\n")
cat(sprintf("   posterior medians from %.4f to %.2f\n", min(med), max(med)))
cat(sprintf("   cows: %d, by number of infected quarters (0/1/2/3/4): %s\n",
            length(med), paste(table(k), collapse = "/")))

# Numbers for the caption
n_cut <- sum(lo < Y_FLOOR)
cat(sprintf("   lower limits cut at %.0e: %d cows; smallest lower limit %.1e\n",
            Y_FLOOR, n_cut, min(lo)))
cat(sprintf("   of which without an infected quarter: %d cows\n", sum(lo < Y_FLOOR & k == 0)))
