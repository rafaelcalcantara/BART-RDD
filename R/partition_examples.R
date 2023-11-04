set.seed(1)
n <- 100
x <- runif(n, -1, 1)
w <- runif(n, -1, 1)
treated <- x >= 0

cex = 0.8
# First split, but opct is not satisfied
png("Figures/bad_trees_1a.png")
plot(x, w, pch = 16 + treated, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
dev.off()

# A split that would not meet the first constraint 
col <- (x > -0.25) & (x < 0) & (w > 0) & (w < 0.4)
png("Figures/bad_trees_1b.png")
plot(x, w, pch = 16 + treated, col = col + 1, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(h = 0.25)
dev.off()

# Splits that are likely to happen to satisfy opct
png("Figures/bad_trees_1c.png")
plot(x, w, pch = 16 + treated, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(v = -0.4)
dev.off()

png("Figures/good_tree_1.png")
plot(x, w, pch = 16 + treated, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(v = -0.4)
abline(v = 0.4)
dev.off()
