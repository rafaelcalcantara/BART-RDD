library(XBART)
library(parallel)
library(foreach)
library(doParallel)
setwd("~/Documents/Git/XBCF-RDD/")
### Function to read results files
readFiles <- function(s,dgp,file)
{
    foreach(x=1:s,.multicombine=T) %do%
        readRDS(paste0("Results/",file,"_dgp",dgp,"_",x,".rds"))
}
#### Parallelization
## no_cores <- detectCores() - 1
## registerDoParallel(no_cores)
## DGP1a
results <- readFiles(s,"1a","cgs")
### Obtain ATE posterior
ate.sum1a <- t(sapply(results, function(x) c(mean(x$fit$atem),quantile(x$fit$atem,c(0.025,0.975)))))
## DGP1b
results <- readFiles(s,"1b","cgs")
### Obtain ATE posterior
ate.sum1b <- t(sapply(results, function(x) c(mean(x$fit$atem),quantile(x$fit$atem,c(0.025,0.975)))))
## DGP2
results <- readFiles(s,"2","cgs")
### Obtain ATE posterior
ate.sum2 <- t(sapply(results, function(x) c(mean(x$fit$atem),quantile(x$fit$atem,c(0.025,0.975)))))
## DGP3
dgp <- readRDS("Data/DGP3.rds")
results <- readFiles(s,"3","cgs")
ate3 <- sapply(dgp,function(x) x$ate) ## Heterogeneous ATE
### Obtain ATE posterior
ate.sum3 <- t(sapply(results, function(x) c(mean(x$fit$atem),quantile(x$fit$atem,c(0.025,0.975)))))
####
res.mat <- cbind(c(mean((ate.sum1a[,1]-0.04)^2),mean(ate.sum1a[,2]<=0.04 & 0.04<=ate.sum1a[,3]),mean(ate.sum1a[,3]-ate.sum1a[,2])),
                 c(mean((ate.sum1b[,1]-0.04)^2),mean(ate.sum1b[,2]<=0.04 & 0.04<=ate.sum1b[,3]),mean(ate.sum1b[,3]-ate.sum1b[,2])),
                 c(mean((ate.sum2[,1]-1)^2),mean(ate.sum2[,2]<=1 & 1<=ate.sum2[,3]),mean(ate.sum2[,3]-ate.sum2[,2])),
                 c(mean((ate.sum3[,1]-ate3)^2),mean(ate.sum3[,2]<=ate3 & ate3<=ate.sum3[,3]),mean(ate.sum3[,3]-ate.sum3[,2])))
colnames(res.mat) <- c("DGP1a","DGP1b","DGP2","DGP3")
rownames(res.mat) <- c("MSE","Coverage","Int. Length")
res.mat
saveRDS(res.mat,"Results/cgs_simulations_1.rds")
write.table(res.mat,"Tables/cgs_simulations_1.csv")
####
## stopImplicitCluster()
####
