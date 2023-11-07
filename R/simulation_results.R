## Reading command line arguments to save results
args <- commandArgs(trailingOnly=TRUE)
sig <- as.numeric(args[1])
k <- as.numeric(args[2])
##
res.all <- readRDS(paste0("Results/results_",sig,"_",k,".rds"))
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
res <- res.all$bart.rdd.1
bart.fit.1 <- res$fit
ate.bart.1 <- colMeans(sapply(bart.fit.1,colMeans))
rmse.bart.1 <- sqrt(mean((ate.bart.1-ate.1)^2))
bias.bart.1 <- mean(ate.bart.1-ate.1)
var.bart.1 <- var(ate.bart.1)
ci.bart.1 <- apply(sapply(bart.fit.1,colMeans),2,quantile,c(0.025,0.975))
ci.bart.1 <- t(ci.bart.1)
cov.bart.1 <- mean(ci.bart.1[,1] <= ate.1 & ate.1 <= ci.bart.1[,2])
size.bart.1 <- mean(ci.bart.1[,2]-ci.bart.1[,1])
#### CATE
ind <- sapply(1:s, function(i) -h.1[i] <= d.1$x[,i] & d.1$x[,i] <= h.1[i])
rmse.cate.bart.1 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bart.fit.1[[i]])-cate.1[ind[,i],i])^2))))
ci.cate.bart.1 <- lapply(1:s,function(i) t(apply(bart.fit.1[[i]],1,quantile,c(0.025,0.975))))
cov.cate.bart.1 <- mean(sapply(1:s, function(i) mean(ci.cate.bart.1[[i]][,1] <= cate.1[ind[,i],i] & cate.1[ind[,i],i] <= ci.cate.bart.1[[i]][,2])))
size.cate.bart.1 <- mean(sapply(1:s, function(i) mean(ci.cate.bart.1[[i]][,2]-ci.cate.bart.1[[i]][,1])))
### DGP2
#### ATE
res <- res.all$bart.rdd.2
bart.fit.2 <- res$fit
ate.bart.2 <- colMeans(sapply(bart.fit.2,colMeans))
rmse.bart.2 <- sqrt(mean((ate.bart.2-ate.2)^2))
bias.bart.2 <- mean(ate.bart.2-ate.2)
var.bart.2 <- var(ate.bart.2)
ci.bart.2 <- apply(sapply(bart.fit.2,colMeans),2,quantile,c(0.025,0.975))
ci.bart.2 <- t(ci.bart.2)
cov.bart.2 <- mean(ci.bart.2[,1] <= ate.2 & ate.2 <= ci.bart.2[,2])
size.bart.2 <- mean(ci.bart.2[,2]-ci.bart.2[,1])
#### CATE
ind <- sapply(1:s, function(i) -h.2[i] <= d.2$x[,i] & d.2$x[,i] <= h.2[i])
rmse.cate.bart.2 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bart.fit.2[[i]])-cate.2[ind[,i],i])^2))))
ci.cate.bart.2 <- lapply(1:s,function(i) t(apply(bart.fit.2[[i]],1,quantile,c(0.025,0.975))))
cov.cate.bart.2 <- mean(sapply(1:s, function(i) mean(ci.cate.bart.2[[i]][,1] <= cate.2[ind[,i],i] & cate.2[ind[,i],i] <= ci.cate.bart.2[[i]][,2])))
size.cate.bart.2 <- mean(sapply(1:s, function(i) mean(ci.cate.bart.2[[i]][,2]-ci.cate.bart.2[[i]][,1])))
### DGP3
#### ATE
res <- res.all$bart.rdd.3
bart.fit.3 <- res$fit
ate.bart.3 <- colMeans(sapply(bart.fit.3,colMeans))
rmse.bart.3 <- sqrt(mean((ate.bart.3-ate.3)^2))
bias.bart.3 <- mean(ate.bart.3-ate.3)
var.bart.3 <- var(ate.bart.3)
ci.bart.3 <- apply(sapply(bart.fit.3,colMeans),2,quantile,c(0.025,0.975))
ci.bart.3 <- t(ci.bart.3)
cov.bart.3 <- mean(ci.bart.3[,1] <= ate.3 & ate.3 <= ci.bart.3[,2])
size.bart.3 <- mean(ci.bart.3[,2]-ci.bart.3[,1])
#### CATE
ind <- sapply(1:s, function(i) -h.3[i] <= d.3$x[,i] & d.3$x[,i] <= h.3[i])
rmse.cate.bart.3 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bart.fit.3[[i]])-cate.3[ind[,i],i])^2))))
ci.cate.bart.3 <- lapply(1:s,function(i) t(apply(bart.fit.3[[i]],1,quantile,c(0.025,0.975))))
cov.cate.bart.3 <- mean(sapply(1:s, function(i) mean(ci.cate.bart.3[[i]][,1] <= cate.3[ind[,i],i] & cate.3[ind[,i],i] <= ci.cate.bart.3[[i]][,2])))
size.cate.bart.3 <- mean(sapply(1:s, function(i) mean(ci.cate.bart.3[[i]][,2]-ci.cate.bart.3[[i]][,1])))
## BCF
### DGP1
#### ATE
res <- res.all$bcf.1
bcf.fit.1 <- res$fit
ate.bcf.1 <- colMeans(sapply(bcf.fit.1,colMeans))
rmse.bcf.1 <- sqrt(mean((ate.bcf.1-ate.1)^2))
bias.bcf.1 <- mean(ate.bcf.1-ate.1)
var.bcf.1 <- var(ate.bcf.1)
ci.bcf.1 <- apply(sapply(bcf.fit.1,colMeans),2,quantile,c(0.025,0.975))
ci.bcf.1 <- t(ci.bcf.1)
cov.bcf.1 <- mean(ci.bcf.1[,1] <= ate.1 & ate.1 <= ci.bcf.1[,2])
size.bcf.1 <- mean(ci.bcf.1[,2]-ci.bcf.1[,1])
#### CATE
rmse.cate.bcf.1 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bcf.fit.1[[i]])-cate.1[,i])^2))))
ci.cate.bcf.1 <- lapply(1:s,function(i) t(apply(bcf.fit.1[[i]],1,quantile,c(0.025,0.975))))
cov.cate.bcf.1 <- mean(sapply(1:s, function(i) mean(ci.cate.bcf.1[[i]][,1] <= cate.1[,i] & cate.1[,i] <= ci.cate.bcf.1[[i]][,2])))
size.cate.bcf.1 <- mean(sapply(1:s, function(i) mean(ci.cate.bcf.1[[i]][,2]-ci.cate.bcf.1[[i]][,1])))
### DGP2
#### ATE
res <- res.all$bcf.2
bcf.fit.2 <- res$fit
ate.bcf.2 <- colMeans(sapply(bcf.fit.2,colMeans))
rmse.bcf.2 <- sqrt(mean((ate.bcf.2-ate.2)^2))
bias.bcf.2 <- mean(ate.bcf.2-ate.2)
var.bcf.2 <- var(ate.bcf.2)
ci.bcf.2 <- apply(sapply(bcf.fit.2,colMeans),2,quantile,c(0.025,0.975))
ci.bcf.2 <- t(ci.bcf.2)
cov.bcf.2 <- mean(ci.bcf.2[,1] <= ate.2 & ate.2 <= ci.bcf.2[,2])
size.bcf.2 <- mean(ci.bcf.2[,2]-ci.bcf.2[,1])
#### CATE
rmse.cate.bcf.2 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bcf.fit.2[[i]])-cate.2[,i])^2))))
ci.cate.bcf.2 <- lapply(1:s,function(i) t(apply(bcf.fit.2[[i]],1,quantile,c(0.025,0.975))))
cov.cate.bcf.2 <- mean(sapply(1:s, function(i) mean(ci.cate.bcf.2[[i]][,1] <= cate.2[,i] & cate.2[,i] <= ci.cate.bcf.2[[i]][,2])))
size.cate.bcf.2 <- mean(sapply(1:s, function(i) mean(ci.cate.bcf.2[[i]][,2]-ci.cate.bcf.2[[i]][,1])))
### DGP3
#### ATE
res <- res.all$bcf.3
bcf.fit.3 <- res$fit
ate.bcf.3 <- colMeans(sapply(bcf.fit.3,colMeans))
rmse.bcf.3 <- sqrt(mean((ate.bcf.3-ate.3)^2))
bias.bcf.3 <- mean(ate.bcf.3-ate.3)
var.bcf.3 <- var(ate.bcf.3)
ci.bcf.3 <- apply(sapply(bcf.fit.3,colMeans),2,quantile,c(0.025,0.975))
ci.bcf.3 <- t(ci.bcf.3)
cov.bcf.3 <- mean(ci.bcf.3[,1] <= ate.3 & ate.3 <= ci.bcf.3[,2])
size.bcf.3 <- mean(ci.bcf.3[,2]-ci.bcf.3[,1])
#### CATE
rmse.cate.bcf.3 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(bcf.fit.3[[i]])-cate.3[,i])^2))))
ci.cate.bcf.3 <- lapply(1:s,function(i) t(apply(bcf.fit.3[[i]],1,quantile,c(0.025,0.975))))
cov.cate.bcf.3 <- mean(sapply(1:s, function(i) mean(ci.cate.bcf.3[[i]][,1] <= cate.3[,i] & cate.3[,i] <= ci.cate.bcf.3[[i]][,2])))
size.cate.bcf.3 <- mean(sapply(1:s, function(i) mean(ci.cate.bcf.3[[i]][,2]-ci.cate.bcf.3[[i]][,1])))
## SBART
### DGP1
#### ATE
res <- res.all$sbart.1
sbart.fit.1 <- res$fit
ate.sbart.1 <- colMeans(sapply(sbart.fit.1,colMeans))
rmse.sbart.1 <- sqrt(mean((ate.sbart.1-ate.1)^2))
bias.sbart.1 <- mean(ate.sbart.1-ate.1)
var.sbart.1 <- var(ate.sbart.1)
ci.sbart.1 <- apply(sapply(sbart.fit.1,colMeans),2,quantile,c(0.025,0.975))
ci.sbart.1 <- t(ci.sbart.1)
cov.sbart.1 <- mean(ci.sbart.1[,1] <= ate.1 & ate.1 <= ci.sbart.1[,2])
size.sbart.1 <- mean(ci.sbart.1[,2]-ci.sbart.1[,1])
#### CATE
rmse.cate.sbart.1 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(sbart.fit.1[[i]])-cate.1[,i])^2))))
ci.cate.sbart.1 <- lapply(1:s,function(i) t(apply(sbart.fit.1[[i]],1,quantile,c(0.025,0.975))))
cov.cate.sbart.1 <- mean(sapply(1:s, function(i) mean(ci.cate.sbart.1[[i]][,1] <= cate.1[,i] & cate.1[,i] <= ci.cate.sbart.1[[i]][,2])))
size.cate.sbart.1 <- mean(sapply(1:s, function(i) mean(ci.cate.sbart.1[[i]][,2]-ci.cate.sbart.1[[i]][,1])))
### DGP2
#### ATE
res <- res.all$sbart.2
sbart.fit.2 <- res$fit
ate.sbart.2 <- colMeans(sapply(sbart.fit.2,colMeans))
rmse.sbart.2 <- sqrt(mean((ate.sbart.2-ate.2)^2))
bias.sbart.2 <- mean(ate.sbart.2-ate.2)
var.sbart.2 <- var(ate.sbart.2)
ci.sbart.2 <- apply(sapply(sbart.fit.2,colMeans),2,quantile,c(0.025,0.975))
ci.sbart.2 <- t(ci.sbart.2)
cov.sbart.2 <- mean(ci.sbart.2[,1] <= ate.2 & ate.2 <= ci.sbart.2[,2])
size.sbart.2 <- mean(ci.sbart.2[,2]-ci.sbart.2[,1])
#### CATE
rmse.cate.sbart.2 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(sbart.fit.2[[i]])-cate.2[,i])^2))))
ci.cate.sbart.2 <- lapply(1:s,function(i) t(apply(sbart.fit.2[[i]],1,quantile,c(0.025,0.975))))
cov.cate.sbart.2 <- mean(sapply(1:s, function(i) mean(ci.cate.sbart.2[[i]][,1] <= cate.2[,i] & cate.2[,i] <= ci.cate.sbart.2[[i]][,2])))
size.cate.sbart.2 <- mean(sapply(1:s, function(i) mean(ci.cate.sbart.2[[i]][,2]-ci.cate.sbart.2[[i]][,1])))
### DGP3
#### ATE
res <- res.all$sbart.3
sbart.fit.3 <- res$fit
ate.sbart.3 <- colMeans(sapply(sbart.fit.3,colMeans))
rmse.sbart.3 <- sqrt(mean((ate.sbart.3-ate.3)^2))
bias.sbart.3 <- mean(ate.sbart.3-ate.3)
var.sbart.3 <- var(ate.sbart.3)
ci.sbart.3 <- apply(sapply(sbart.fit.3,colMeans),2,quantile,c(0.025,0.975))
ci.sbart.3 <- t(ci.sbart.3)
cov.sbart.3 <- mean(ci.sbart.3[,1] <= ate.3 & ate.3 <= ci.sbart.3[,2])
size.sbart.3 <- mean(ci.sbart.3[,2]-ci.sbart.3[,1])
#### CATE
rmse.cate.sbart.3 <- mean(sapply(1:s,function(i) sqrt(mean((rowMeans(sbart.fit.3[[i]])-cate.3[,i])^2))))
ci.cate.sbart.3 <- lapply(1:s,function(i) t(apply(sbart.fit.3[[i]],1,quantile,c(0.025,0.975))))
cov.cate.sbart.3 <- mean(sapply(1:s, function(i) mean(ci.cate.sbart.3[[i]][,1] <= cate.3[,i] & cate.3[,i] <= ci.cate.sbart.3[[i]][,2])))
size.cate.sbart.3 <- mean(sapply(1:s, function(i) mean(ci.cate.sbart.3[[i]][,2]-ci.cate.sbart.3[[i]][,1])))
## TBART
### DGP1
#### ATE
res <- res.all$tbart.1
tbart.fit.1 <- res$fit
ate.tbart.1 <- sapply(tbart.fit.1,mean)
rmse.tbart.1 <- sqrt(mean((ate.tbart.1-ate.1)^2))
bias.tbart.1 <- mean(ate.tbart.1-ate.1)
var.tbart.1 <- var(ate.tbart.1)
ci.tbart.1 <- sapply(tbart.fit.1,quantile,c(0.025,0.975))
ci.tbart.1 <- t(ci.tbart.1)
cov.tbart.1 <- mean(ci.tbart.1[,1] <= ate.1 & ate.1 <= ci.tbart.1[,2])
size.tbart.1 <- mean(ci.tbart.1[,2]-ci.tbart.1[,1])
### DGP2
#### ATE
res <- res.all$tbart.2
tbart.fit.2 <- res$fit
ate.tbart.2 <- sapply(tbart.fit.2,mean)
rmse.tbart.2 <- sqrt(mean((ate.tbart.2-ate.2)^2))
bias.tbart.2 <- mean(ate.tbart.2-ate.2)
var.tbart.2 <- var(ate.tbart.2)
ci.tbart.2 <- sapply(tbart.fit.2,quantile,c(0.025,0.975))
ci.tbart.2 <- t(ci.tbart.2)
cov.tbart.2 <- mean(ci.tbart.2[,1] <= ate.2 & ate.2 <= ci.tbart.2[,2])
size.tbart.2 <- mean(ci.tbart.2[,2]-ci.tbart.2[,1])
### DGP3
#### ATE
res <- res.all$tbart.3
tbart.fit.3 <- res$fit
ate.tbart.3 <- sapply(tbart.fit.3,mean)
rmse.tbart.3 <- sqrt(mean((ate.tbart.3-ate.3)^2))
bias.tbart.3 <- mean(ate.tbart.3-ate.3)
var.tbart.3 <- var(ate.tbart.3)
ci.tbart.3 <- sapply(tbart.fit.3,quantile,c(0.025,0.975))
ci.tbart.3 <- t(ci.tbart.3)
cov.tbart.3 <- mean(ci.tbart.3[,1] <= ate.3 & ate.3 <= ci.tbart.3[,2])
size.tbart.3 <- mean(ci.tbart.3[,2]-ci.tbart.3[,1])
## CKT
### DGP1
res <- res.all$ckt.1
ckt.fit.1 <- res$fit
ate.ckt.1 <- sapply(ckt.fit.1, function(i) i$coef[3])
rmse.ckt.1 <- sqrt(mean((ate.ckt.1-ate.1)^2))
bias.ckt.1 <- mean(ate.ckt.1-ate.1)
var.ckt.1 <- var(ate.ckt.1)
ci.ckt.1 <- sapply(ckt.fit.1, function(i) i$ci[3,])
ci.ckt.1 <- t(ci.ckt.1)
cov.ckt.1 <- mean(ci.ckt.1[,1] <= ate.1 & ate.1 <= ci.ckt.1[,2])
size.ckt.1 <- mean(ci.ckt.1[,2]-ci.ckt.1[,1])
### DGP2
res <- res.all$ckt.2
ckt.fit.2 <- res$fit
ate.ckt.2 <- sapply(ckt.fit.2, function(i) i$coef[3])
rmse.ckt.2 <- sqrt(mean((ate.ckt.2-ate.2)^2))
bias.ckt.2 <- mean(ate.ckt.2-ate.2)
var.ckt.2 <- var(ate.ckt.2)
ci.ckt.2 <- sapply(ckt.fit.1, function(i) i$ci[3,])
ci.ckt.2 <- t(ci.ckt.2)
cov.ckt.2 <- mean(ci.ckt.2[,1] <= ate.2 & ate.2 <= ci.ckt.2[,2])
size.ckt.2 <- mean(ci.ckt.2[,2]-ci.ckt.2[,1])
### DGP3
res <- res.all$ckt.3
ckt.fit.3 <- res$fit
ate.ckt.3 <- sapply(ckt.fit.3, function(i) i$coef[3])
rmse.ckt.3 <- sqrt(mean((ate.ckt.3-ate.3)^2))
bias.ckt.3 <- mean(ate.ckt.3-ate.3)
var.ckt.3 <- var(ate.ckt.3)
ci.ckt.3 <- sapply(ckt.fit.3, function(i) i$ci[3,])
ci.ckt.3 <- t(ci.ckt.3)
cov.ckt.3 <- mean(ci.ckt.3[,1] <= ate.3 & ate.3 <= ci.ckt.3[,2])
size.ckt.3 <- mean(ci.ckt.3[,2]-ci.ckt.3[,1])
## Compiling results
ate.res <- data.frame(BART.RDD=c(rmse.bart.1,bias.bart.1,var.bart.1,cov.bart.1,size.bart.1,
                                 rmse.bart.2,bias.bart.2,var.bart.2,cov.bart.2,size.bart.2,
                                 rmse.bart.3,bias.bart.3,var.bart.3,cov.bart.3,size.bart.3),
                      BCF=c(rmse.bcf.1,bias.bcf.1,var.bcf.1,cov.bcf.1,size.bcf.1,
                            rmse.bcf.2,bias.bcf.2,var.bcf.2,cov.bcf.2,size.bcf.2,
                            rmse.bcf.3,bias.bcf.3,var.bcf.3,cov.bcf.3,size.bcf.3),
                      SBART=c(rmse.sbart.1,bias.sbart.1,var.sbart.1,cov.sbart.1,size.sbart.1,
                              rmse.sbart.2,bias.sbart.2,var.sbart.2,cov.sbart.2,size.sbart.2,
                              rmse.sbart.3,bias.sbart.3,var.sbart.3,cov.sbart.3,size.sbart.3),
                      TBART=c(rmse.tbart.1,bias.tbart.1,var.tbart.1,cov.tbart.1,size.tbart.1,
                              rmse.tbart.2,bias.tbart.2,var.tbart.2,cov.tbart.2,size.tbart.2,
                              rmse.tbart.3,bias.tbart.3,var.tbart.3,cov.tbart.3,size.tbart.3),
                      CKT=c(rmse.ckt.1,bias.ckt.1,var.ckt.1,cov.ckt.1,size.ckt.1,
                            rmse.ckt.2,bias.ckt.2,var.ckt.2,cov.ckt.2,size.ckt.2,
                            rmse.ckt.3,bias.ckt.3,var.ckt.3,cov.ckt.3,size.ckt.3))
rownames(ate.res) <- paste(c("RMSE","Bias","Variance","Coverage","Size"),c(rep(1,5),rep(2,5),rep(3,5)),sep=".")
saveRDS(round(ate.res,3),paste0("Tables/ate_res_",sig,"_",k,".rds"))
###
cate.res <- data.frame(BART.RDD=c(rmse.cate.bart.1,cov.cate.bart.1,size.cate.bart.1,
                                  rmse.cate.bart.2,cov.cate.bart.2,size.cate.bart.2,
                                  rmse.cate.bart.3,cov.cate.bart.3,size.cate.bart.3),
                       BCF=c(rmse.cate.bcf.1,cov.cate.bcf.1,size.cate.bcf.1,
                             rmse.cate.bcf.2,cov.cate.bcf.2,size.cate.bcf.2,
                             rmse.cate.bcf.3,cov.cate.bcf.3,size.cate.bcf.3),
                       SBART=c(rmse.cate.sbart.1,cov.cate.sbart.1,size.cate.sbart.1,
                               rmse.cate.sbart.2,cov.cate.sbart.2,size.cate.sbart.2,
                               rmse.cate.sbart.3,cov.cate.sbart.3,size.cate.sbart.3))
rownames(cate.res) <- paste(c("RMSE","Coverage","Size"),c(rep(1,3),rep(2,3),rep(3,3)),sep=".")
saveRDS(round(cate.res,3),paste0("Tables/cate_res_",sig,"_",k,".rds"))
##
print("Done!")
