set.seed(1)
setwd("~/../Git/BART-RDD/")
n <- 100
x <- runif(n, -1, 1)
w <- runif(n, -1, 1)
treated <- x >= 0

cex = 0.8
# First split, but opct is not satisfied
pdf("Figures/bad_trees_1a.pdf")
plot(x, w, pch = 21 + 3*treated,bg=ifelse(treated,"black","white"),col=1, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
dev.off()

# A split that would not meet the first constraint 
col <- (x > -0.25) & (x < 0) & (w > 0) & (w < 0.4)
pdf("Figures/bad_trees_1b.pdf")
plot(x, w, pch = 21 + 3*treated,bg=ifelse(treated,"black",ifelse(col,"red","white")),col=1, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(h = 0.25)
dev.off()

# Splits that are likely to happen to satisfy opct
pdf("Figures/bad_trees_1c.pdf")
plot(x, w, pch = 21 + 3*treated,bg=ifelse(treated,"black","white"),col=1, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(v = -0.4)
dev.off()

pdf("Figures/good_tree_1.pdf")
plot(x, w, pch = 21 + 3*treated,bg=ifelse(treated,"black","white"),col=1, cex = cex)
abline(v = -0.25, lty = 2)
abline(v = 0, lty = 2)
abline(v = 0.25, lty = 2)
abline(h = 0)
abline(v = -0.4)
abline(v = 0.4)
dev.off()
