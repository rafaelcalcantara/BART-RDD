## Data
d.1 <- readRDS("Data/dgp1.rds")
d.2 <- readRDS("Data/dgp2.rds")
d.3 <- readRDS("Data/dgp3.rds")
h.1 <- apply(d.1$x,2,sd)*0.5
h.2 <- apply(d.2$x,2,sd)*0.5
h.3 <- apply(d.3$x,2,sd)*0.5
s <- ncol(d.1$x)
cate.1 <- d.1$cate
ate.1 <- colMeans(cate.1)
cate.2 <- d.2$cate
ate.2 <- colMeans(cate.2)
cate.3 <- d.3$cate
ate.3 <- colMeans(cate.3)
## BART-RDD
### DGP1
#### ATE
res <- readRDS("Results/bart_rdd_1.rds")
bart.fit.1 <- res$fit
ate.bart.1 <- colMeans(sapply(bart.fit.1,colMeans))
rmse.bart.1 <- sqrt(mean((ate.bart.1-ate.1)^2))
bias.bart.1 <- mean(ate.bart.1-ate.1)
var.bart.1 <- var(ate.bart.1)
ci.bart.1 <- apply(sapply(bart.fit.1,colMeans),2,quantile,c(0.025,0.975))
ci.bart.1 <- t(ci.bart.1)
cov.bart.1 <- mean(sapply(1:s,function(i) ci.bart.1[,1] <= ate.1[i] & ate.1[i] <= ci.bart.1[,2]))
#### CATE
ind <- sapply(1:s, function(i) -h.1[i] <= d.1$x[,i] & d.1$x[,i] <= h.1[i])
rmse.cate.1 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bart.fit.1[[i]])-cate.1[ind[,i]])^2))))
ci.cate.1 <- lapply(1:s,function(i) t(apply(bart.fit.1[[i]],1,quantile,c(0.025,0.975))))
cov.cate.1 <- mean(sapply(1:s, function(i) mean(ci.cate.1[[i]][,1] <= cate.1[ind[,i],i] & cate.1[ind[,i],i] <= ci.cate.1[[i]][,2])))
### DGP2
#### ATE
res <- readRDS("Results/bart_rdd_2.rds")
bart.fit.2 <- res$fit
ate.bart.2 <- colMeans(sapply(bart.fit.2,colMeans))
rmse.bart.2 <- sqrt(mean((ate.bart.2-ate.2)^2))
bias.bart.2 <- mean(ate.bart.2-ate.2)
var.bart.2 <- var(ate.bart.2)
ci.bart.2 <- apply(sapply(bart.fit.2,colMeans),2,quantile,c(0.025,0.975))
ci.bart.2 <- t(ci.bart.2)
cov.bart.2 <- mean(sapply(1:s,function(i) ci.bart.2[,1] <= ate.2[i] & ate.2[i] <= ci.bart.2[,2]))
#### CATE
ind <- sapply(1:s, function(i) -h.2[i] <= d.2$x[,i] & d.2$x[,i] <= h.2[i])
rmse.cate.2 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bart.fit.2[[i]])-cate.2[ind[,i]])^2))))
ci.cate.2 <- lapply(1:s,function(i) t(apply(bart.fit.2[[i]],1,quantile,c(0.025,0.975))))
cov.cate.2 <- mean(sapply(1:s, function(i) mean(ci.cate.2[[i]][,1] <= cate.2[ind[,i],i] & cate.2[ind[,i],i] <= ci.cate.2[[i]][,2])))
### DGP3
#### ATE
res <- readRDS("Results/bart_rdd_3.rds")
bart.fit.3 <- res$fit
ate.bart.3 <- colMeans(sapply(bart.fit.3,colMeans))
rmse.bart.3 <- sqrt(mean((ate.bart.3-ate.3)^2))
bias.bart.3 <- mean(ate.bart.3-ate.3)
var.bart.3 <- var(ate.bart.3)
ci.bart.3 <- apply(sapply(bart.fit.3,colMeans),2,quantile,c(0.025,0.975))
ci.bart.3 <- t(ci.bart.3)
cov.bart.3 <- mean(sapply(1:s,function(i) ci.bart.3[,1] <= ate.3[i] & ate.3[i] <= ci.bart.3[,2]))
#### CATE
ind <- sapply(1:s, function(i) -h.3[i] <= d.3$x[,i] & d.3$x[,i] <= h.3[i])
rmse.cate.3 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bart.fit.3[[i]])-cate.3[ind[,i]])^2))))
ci.cate.3 <- lapply(1:s,function(i) t(apply(bart.fit.3[[i]],1,quantile,c(0.025,0.975))))
cov.cate.3 <- mean(sapply(1:s, function(i) mean(ci.cate.3[[i]][,1] <= cate.3[ind[,i],i] & cate.3[ind[,i],i] <= ci.cate.3[[i]][,2])))
