## Join results for each method
setwd("~/Documents/Git/XBCF-RDD")
xbcf <- readRDS("Results/results_xbcf.rds")
cgs <- readRDS("Results/results_cgs.rds")
kr <- readRDS("Results/results_kr.rds")
results <- list(`XBCF-RDD`=xbcf,CGS=cgs,KR=kr)
## MSE
mse <- do.call("rbind",lapply(results, function(x) do.call("cbind",lapply(x, function(y) mean((y[,1]-y[,4])^2)))))
rownames(mse) <- c("XBCF-RDD","CGS","KR")
colnames(mse) <- c("DGP 1 (a)", "DGP 1 (b)", "DGP 2", "DGP 3")
## Variance
Var <- do.call("rbind",lapply(results, function(x) do.call("cbind",lapply(x, function(y) var(y[,1])))))
rownames(Var) <- c("XBCF-RDD","CGS","KR")
colnames(Var) <- c("DGP 1 (a)", "DGP 1 (b)", "DGP 2", "DGP 3")
## Interval length
int.length <- do.call("rbind",lapply(results, function(x) do.call("cbind",lapply(x, function(y) mean(y[,3]-y[,2])))))
rownames(int.length) <- c("XBCF-RDD","CGS","KR")
colnames(int.length) <- c("DGP 1 (a)", "DGP 1 (b)", "DGP 2", "DGP 3")
## Coverage rate
cov.rate <- do.call("rbind",lapply(results, function(x) do.call("cbind",lapply(x, function(y) mean(y[,4] >= y[,2] & y[,4] <= y[,3])))))
rownames(cov.rate) <- c("XBCF-RDD","CGS","KR")
colnames(cov.rate) <- c("DGP 1 (a)", "DGP 1 (b)", "DGP 2", "DGP 3")
