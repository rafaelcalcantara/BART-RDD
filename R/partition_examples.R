set.seed(1)
n <- 100
x <- runif(n, -1, 1)
w <- runif(n, -1, 1)
treated <- x >= 0

par(mfrow = c(2, 2))
cex = 0.8
# First split, but opct is not satisfied
plot(x, w, pch = 16 + treated, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)

# A split that would not meat the first constraint 
col <- (x > -0.25) & (x < 0) & (w > 0) & (w < 0.4)
plot(x, w, pch = 16 + treated, col = col + 1, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(h = 0.25)

# Splits that are likely to happen to satisfy opct
plot(x, w, pch = 16 + treated, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(v = -0.4)

plot(x, w, pch = 16 + treated, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(v = -0.4)
abline(v = 0.4)
