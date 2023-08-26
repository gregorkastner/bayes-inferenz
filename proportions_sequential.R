# Beta-Binomial Model to illustrate some things in the book chapter

set.seed(42)

# Illustrate sequential updating:

a0 <- 1
b0 <- 1
N <- 100 
thetas <- c(0, .1, .5)

xs <- seq(0, 1, length.out = 201)

dir.create('plots', showWarnings = FALSE)

pdf('plots/prior2posterior.pdf', width = 8, height = 3.5)
par(mar = c(1.6, 1.6, .5, .5), mgp = c(1.6, .6, 0), mfrow = c(1,3))
for (th in seq_along(thetas)) {
  plot(xs, dbeta(xs, a0, b0), type = 'l', ylim = c(0, 11), col = rgb(0,0,0,.2),
       xlab = '', ylab = '', main = '')
  legend('topright', paste("N =", c(0, 20, 40, 60, 80, 100)), lty = 1,
         col = rgb(0, 0, 0, .2 + .4*c(0, 20, 40, 60, 80, 100)/N))
  successes <- failures <- 0L
  for (i in seq_len(N)) {
    if (rbinom(1, 1, thetas[th])) successes <- successes + 1L else failures <- failures + 1L
    lines(xs, dbeta(xs, a0 + successes, b0 + failures), col = rgb(0, 0, 0, .2 + .4*i/N))
  }
}
dev.off()
