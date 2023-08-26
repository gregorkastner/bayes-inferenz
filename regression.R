# Regression illustration

set.seed(42)

library(stochvol) # contains exchange rate data
library(mvtnorm) # for dmvt and rmvnorm
library(rgl) # for posterior visualisation
library(invgamma) # for quantiles of the inverse gamma distribution

data(exrates)

dat <- log(exrates$USD)
dates <- exrates$date

dir.create('plots', showWarnings = FALSE)

pdf('plots/regData.pdf', width = 10, height = 4)
par(mar = c(1.5, 1.5, .5, .5), mgp = c(1.6, .6, 0))
plot(exp(dat), xaxt = 'n', type = 'l', xlab = '', ylab = '', log = 'y')
ats = seq(1, length(dates), length.out = 8)
axis(1, at = ats, labels = dates[ats])
dev.off()

# Construct the design matrix X and the response vector y
X <- cbind(1, head(dat, -1))
y <- tail(dat, -1)
N <- length(y)

# Compute posterior statistics
betahat <- tcrossprod(solve(crossprod(X)), X) %*% y
errors <- y - X %*% betahat
SSR <- sum(errors^2)
df <- nrow(X) - ncol(X)
bN <- betahat
BN <- solve(crossprod(X))
cN <- df/2
CN <- SSR/2
sigma2hat <- SSR/df

# Generate posterior draws
ndraws <- 100000
sigma2draws <- 1/rgamma(ndraws, cN, CN)
sigmadraws <- sqrt(sigma2draws)
betadraws <- t(as.numeric(bN) + t(sigmadraws * rmvnorm(ndraws, rep(0, ncol(X)), BN)))
b1draws <- betadraws[,1]
b2draws <- betadraws[,2]

# Visualize
gridlen <- 101
b1grid <- seq(quantile(b1draws, .05), quantile(b1draws, .95), length.out = gridlen)
b2grid <- seq(quantile(b2draws, .05), quantile(b2draws, .95), length.out = gridlen)
bgrid <- expand.grid(b1grid, b2grid)
density = matrix(dmvt(bgrid, delta = bN, sigma = sigma2hat * BN, df = 2*cN, log = FALSE), nrow = gridlen)
#density = matrix(dmvnorm(bgrid, bN, sigma2hat * BN, log = FALSE), nrow = gridlen)
#persp3d(b1grid, b2grid, density, col = heat.colors(length(density))[rank(density)])
#rgl.snapshot('plots/regression_posterior_bivariate.png')


pdf('plots/regPostBeta.pdf', width = 10, height = 4)
par(mar = c(2, 4, .5, .5), mgp = c(1.6, .6, 0))
filled.contour(b1grid, b2grid, density, plot.axes = {
  axis(1)
  axis(2)
  contour(b1grid, b2grid, density, add = TRUE, lwd = 2)
  })
dev.off()

xgrid <- seq(-.03, .02, length.out = gridlen)
Xgrid <- cbind(1, xgrid)
ygrid <- mean(b1draws) + mean(b2draws) * xgrid
xygrid <- expand.grid(xgrid, ygrid)
condmeans <- as.numeric(Xgrid %*% bN)
condvars <- (rowSums(Xgrid %*% BN * Xgrid) + 1) * CN / cN

densit1 <- matrix(dt((xygrid[,2] - condmeans)/sqrt(condvars), df = 2*cN) / sqrt(condvars), nrow = gridlen, byrow = TRUE)

pdf('plots/regPredXY.pdf', width = 10, height = 4)
par(mar = c(2.5, 3.5, .5, .5), mgp = c(1.6, .6, 0))

filled.contour(xgrid, ygrid, densit1, plot.axes = {
  axis(1)
  mtext(expression(y[N]), 1, 1.3)
  axis(2)
  mtext(expression(y[N+1]), 2, 2)
  points(X[,2], y)
  })

dev.off()

# Credible intervals:
alpha <- .05
beta1CI <- bN[1] + sqrt(sigma2hat * BN[1,1]) * qt(c(alpha/2, 1 - alpha/2), df = 2*cN)
beta2CI <- bN[2] + sqrt(sigma2hat * BN[2,2]) * qt(c(alpha/2, 1 - alpha/2), df = 2*cN)
sigma2CI <- qinvgamma(c(alpha/2, 1 - alpha/2),  cN, CN)

# Compute HPD intervals (just to check)
resolution <- 100000
grid <- seq(0, 1, length.out = resolution + 1)
dist <- (1 - alpha) * resolution

# beta1
qs <- bN[1] + sqrt(sigma2hat * BN[1,1]) * qt(grid, df = 2*cN)
minimizer <- which.min(diff(qs, lag = dist))
beta1HPD <- c(qs[minimizer], qs[minimizer + dist])

print(pt((0 - bN[1])/sqrt(sigma2hat * BN[1,1]), df = 2*cN))
print(pt((1 - bN[2])/sqrt(sigma2hat * BN[2,2]), df = 2*cN, lower.tail = FALSE))

# beta2
qs <- bN[2] + sqrt(sigma2hat * BN[2,2]) * qt(grid, df = 2*cN)
minimizer <- which.min(diff(qs, lag = dist))
beta2HPD <- c(qs[minimizer], qs[minimizer + dist])

# sigma2
qs <- qinvgamma(grid, cN, CN)
minimizer <- which.min(diff(qs, lag = dist))
sigma2HPD <- c(qs[minimizer], qs[minimizer + dist])

pdf('plots/regPostSigma2.pdf', width = 5, height = 5)
s2 <- seq(quantile(sigma2draws, 0.001), quantile(sigma2draws, .999), length.out = gridlen)
par(mar = c(1.5, 1.5, 1.5, .5), mgp = c(1.6, .6, 0))
plot(s2, dinvgamma(s2, cN, CN), type = 'l', xlab = '', ylab = '')
dev.off()

# Compute the long run mean (after restricting the prior to enforce stationarity)
pdf('plots/regPostMean.pdf', width = 12, height = 4)
par(mar = c(1.5, 1.5, 0, 0), mgp = c(1.6, .6, 0), mfrow = c(1, 3))
selecta <- abs(b2draws) <= .9999
hist(b2draws[selecta], breaks = 42, main = '', prob = TRUE)
hist(b1draws[selecta] / (1 - b2draws[selecta]), breaks = 42, main = '', prob = TRUE)
hist(sqrt(sigma2draws[selecta] / (1 - b2draws[selecta]^2)), breaks = 42, main = '', prob = TRUE)
dev.off()

# Multi-Step predictions:
len <- sum(selecta)
ahead <- 5
meandraws <- preddraws <- matrix(NA_real_, len, ahead)
sddraws <- sqrt(sigma2draws[selecta])

meandraws[,1] <- b1draws[selecta] + b2draws[selecta] * tail(y, 1)
preddraws[,1] <- rnorm(len, meandraws[,1], sddraws)
for (i in 2:5) {
  meandraws[,i] <- b1draws[selecta] + b2draws[selecta] * preddraws[,i-1]
  preddraws[,i] <- rnorm(len, meandraws[,i], sddraws)
}

adjust <- 2
pdf('plots/ARpred.pdf', width = 10, height = 4)
par(mar = c(1.6, 1.6, .5, .5), mgp = c(1.6, .6, 0))
plot(density(preddraws[,1], bw = "SJ", adjust = adjust),
     main = '', xlab = '', ylab = '', lwd = 2)
for (i in 2:5) {
  lines(density(preddraws[,i], bw = "SJ", adjust = adjust), col = i, lwd = 2, lty = i)
}
legend("topright", paste("h = ", 1:5), col = 1:ahead, lwd = 2, lty = 1:5)
dev.off()

# Find the maxima:
whichmaxs <- apply(preddraws, 1, which.max)
maxs <- apply(preddraws, 1, max)

pdf('plots/ARmax.pdf', width = 10, height = 4)
par(mar = c(1.6, 1.6, .5, .5), mgp = c(1.6, .6, 0), mfrow = c(1, 2))
hist(maxs, breaks = 40, xlab = '', ylab = '', main = '', probability = TRUE)
barplot(table(whichmaxs)/length(whichmaxs), main = '')
dev.off()


# Testing:
# Prior hyperparameters:
c0 <- 11
C0 <- 1
C0/(c0-1)
sqrt(C0^2/((c0-1)^2 * (c0-2)))
b0 <- c(0,1)
B0 <- 1*diag(2)

Xorig <- X
b0orig <- b0
B0orig <- B0
yorig <- y
dropintercept <- FALSE
dropAR <- FALSE

# log ML of model without any predictors
logML0 <- (-N/2) * log(2*pi) + lgamma(N/2 + c0) - lgamma(c0) + c0 * log(C0) - (N/2 + c0) * log(sum(diff(dat)^2)/2 + C0)

for (i in 1:3) {
  if (i == 2) dropintercept <- TRUE
  else if (i == 3) {
    dropintercept <- FALSE
    dropAR <- TRUE
  }

if (dropintercept) {
  X <- Xorig[,-1, drop = FALSE]
  b0 <- b0orig[2]
  B0 <- matrix(B0orig[2,2])
}

if (dropAR) {
  X <- Xorig[,-2, drop = FALSE]
  y <- diff(dat)
  b0 <- b0orig[1]
  B0 <- matrix(B0orig[1,1])
}

# Compute posterior statistics
B0inv <- solve(B0)
BNinv <- B0inv + crossprod(X)
BN <- solve(BNinv)
bN <- BN %*% (B0inv %*% b0 + crossprod(X, y))
err <- as.numeric(crossprod(y) + crossprod(b0, B0inv %*% b0) - crossprod(bN, BNinv %*% bN))
cN <- c0 + N/2
CN <- C0 + err/2
sigma2Ihat <- CN/cN

# log ML:
if (dropintercept) {
  logML2 <- (-N/2) * log(2*pi) + .5 * as.numeric(determinant(B0inv, log = TRUE)$modulus) - as.numeric(determinant(BNinv, log = TRUE)$modulus) + c0*log(C0) - cN *log(CN) + lgamma(cN) - lgamma(c0)  
  } else if (dropAR) {
  logML1 <- (-N/2) * log(2*pi) + .5 * as.numeric(determinant(B0inv, log = TRUE)$modulus) - as.numeric(determinant(BNinv, log = TRUE)$modulus) + c0*log(C0) - cN *log(CN) + lgamma(cN) - lgamma(c0)
} else {
  logML12 <- (-N/2) * log(2*pi) + .5 * as.numeric(determinant(B0inv, log = TRUE)$modulus) - as.numeric(determinant(BNinv, log = TRUE)$modulus) + c0*log(C0) - cN *log(CN) + lgamma(cN) - lgamma(c0)
  # Savage-Dickey is biased:
  logSD0v12 <- dmvt(c(0,1), delta = b0, sigma = C0/c0 * B0, df = 2*c0, log = TRUE) - dmvt(c(0,1), delta = bN, sigma = sigma2Ihat * BN, df = 2*cN, log = TRUE)
}
}

logMLS <- c(logML0, logML1, logML2, logML12)
logMLS - logMLS[1]
logSD0v12
