# Beta-Binomial Model to illustrate some things in the book chapter

set.seed(42)

# Set simdata to FALSE to use the data from the first SORA prevalence study:
# https://www.sora.at/fileadmin/downloads/projekte/Austria_Spread_of_SARS-CoV-2_Study_Report.pdf

simdata <- FALSE

if (simdata) {
  # DGP:
  theta_true <- .1
  n <- 100
  # Generate the data:
  data <- sample(c(0, 1), n, replace = TRUE, prob = c(1 - theta_true, theta_true))
  successes <- sum(data)
  failures <- n - successes
} else {
  n <- 1544
  successes <- 6
  failures <- n - successes
}

# Fixed values:
theta0_1 <- .005
theta0_2 <- .01

# Priors:
prior_mean <- c(.5, .5, .01)
prior_sd <- c(1/sqrt(12), 1/sqrt(8), .01)

prior_a <- ((1 - prior_mean) / prior_sd^2 - 1 / prior_mean) * prior_mean^2
prior_b <- prior_a * (1 / prior_mean - 1)

# Posteriors:
post_a <- successes + prior_a
post_b <- failures + prior_b

dir.create('plots', showWarnings = FALSE)

# Visualization of priors and posteriors:
x <- seq(0.0001, .02, length.out = 200)
pdf('plots/prop.pdf', width = 10, height = 4)
par(mgp = c(1.6, .6, 0), mar = c(2, 2, .5, .5), mfrow = c(1, 2))
for (logscale in c('', 'y')) {
  plot(x, dbeta(x, post_a[1], post_b[1]), col = 0, xlab = '', ylab = '', log = logscale)
  for (i in 1:3) {
    lines(x, dbeta(x, post_a[i], post_b[i]), type = 'l', lwd = 2, xlab = '', ylab = '', cex = 2, col = i)
    lines(x, dbeta(x, prior_a[i], prior_b[i]), col = i, lty = 2, lwd = 2, pch = 2, cex = 2)
  }
}
dev.off()

# "Symmetrical" 95% credible intervals
alpha <- .05
left <- qbeta(alpha/2, post_a, post_b)
right <- qbeta(1 - alpha/2, post_a, post_b)

# HPDs
resolution <- 10000
grid <- seq(0, 1, length.out = resolution + 1)
dist <- (1 - alpha) * resolution

leftHPD <- rightHPD <- rep(NA_real_, 3)
for (i in 1:3) {
    qs <- qbeta(grid, post_a[i], post_b[i])
    minimizer <- which.min(diff(qs, lag = dist))
    leftHPD[i] <- qs[minimizer]
    rightHPD[i] <- qs[minimizer + dist]
}

# Posterior Predictive:
beta_binom <- function(k, H, aN, bN, log = FALSE) {
  lval <- lchoose(H, k) + lbeta(aN + k, bN + H - k) - lbeta(aN, bN)
  if (log) lval else exp(lval)
}

predvals <- 0:5
H <- 100

pdf('plots/predprop.pdf', width = 10, height = 4)
par(mgp = c(1.6, .6, 0), mar = c(2, 2, .5, .5), mfrow = c(1,2))
plot(predvals, beta_binom(predvals, H, post_a[3], post_b[3]), type = 'b', lwd = 2, xlab = '', ylab = '', cex = 2)
lines(predvals, dbinom(predvals, H, post_a[3] / (post_a[3] + post_b[3])), type = 'b',
      col = 2, lty = 2, lwd = 2, pch = 2, cex = 2)
plot(predvals, beta_binom(predvals, H, post_a[3], post_b[3]), type = 'b', lwd = 2, xlab = '', ylab = '', cex = 2, log = 'y', ylim = range(dbinom(predvals, H, post_a[3] / (post_a[3] + post_b[3]))))
lines(predvals, dbinom(predvals, H, post_a[3] / (post_a[3] + post_b[3])), type = 'b',
      col = 2, lty = 2, lwd = 2, pch = 2, cex = 2)
dev.off()

# Posterior densitity:
post_dens_unnormalized <- function(theta, n, successes, a, b, log = FALSE) {
  tmp <- dbinom(successes, size = n, prob = theta, log = TRUE) +
    dbeta(theta, shape1 = a, shape2 = b, log = TRUE)
  if(log) tmp else exp(tmp)
}

# Log marginal likelihoods:
lmarglik_fixed <- dbinom(successes, size = n, prob = c(theta0_1, theta0_2), log = TRUE)

lmarglik <- log(integrate(post_dens_unnormalized, lower = 0, upper = 1, n = n,
  successes = successes, a = prior_a[3], b = prior_b[3])$value)

# analytically:
lmarglik <- lchoose(n, successes) - lbeta(prior_a[3], prior_b[3]) +
    lbeta(successes + prior_a[3], failures + prior_b[3])

# Log Bayes factors
print(lBF <- lmarglik_fixed - lmarglik)

# Log Savage-Dickey density ratios
lSD <- dbeta(c(theta0_1, theta0_2), post_a[3], post_b[3], log = TRUE) - 
  dbeta(c(theta0_1, theta0_2), prior_a[3], prior_b[3], log = TRUE)

# Illustration:
xgrid <- seq(0, .02, length.out = 200)
lwd <- 2
pdf('plots/savage-dickey.pdf', width = 12, height = 5)
par(mgp = c(1.6, .6, 0), mar = c(2, 2, .5, .5))
plot(xgrid, dbeta(xgrid, post_a[3], post_b[3]), type = 'l', xlab = '', ylab = '', lwd = lwd)
lines(xgrid, dbeta(xgrid, prior_a[3], prior_b[3]), col = 3, lty = 2, lwd = 1.5 * lwd)

for (theta0 in c(theta0_1, theta0_2)) {
  lines(rep(theta0, 2), c(0, dbeta(theta0, post_a[3], post_b[3])), lwd = lwd)
  lines(rep(theta0, 2), c(0, dbeta(theta0, prior_a[3], prior_b[3])), col = 3, lty = 2, lwd = 1.5 * lwd)
  points(theta0, 0, col = 2, cex = 3, lwd = lwd)
}

abline(h = 0, lwd = 1)
dev.off()
