library(XBART)
library(parallel)
library(foreach)
library(doParallel)
setwd("~/Documents/Git/XBCF-RDD/")
### Function to read results files
readFiles <- function(s,dgp,file)
{
    foreach(x=1:s,.multicombine=T) %dopar%
        readRDS(paste0("Results/",file,"_dgp",dgp,"_",x,".rds"))
}
#### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## DGP1a
results <- readFiles(s,"1a","xbcf")
### Obtain ATE posterior
ate.post1a <- t(sapply(results,function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
## DGP1b
results <- readFiles(s,"1b","xbcf")
### Obtain ATE posterior
ate.post1b <- t(sapply(results,function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
## DGP2
results <- readFiles(s,"2","xbcf")
### Obtain ATE posterior
ate.post2 <- t(sapply(results,function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
## DGP3
dgp <- readRDS("Data/DGP3.rds")
results <- readFiles(s,"3","xbcf")
ate3 <- sapply(dgp,function(x) x$ate) ## Heterogeneous ATE
### Obtain ATE posterior
ate.post3 <- t(sapply(results,function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
####
res.mat <- cbind(c(mean((ate.post1a[,1]-0.04)^2),mean(ate.post1a[,2]<=0.04 & 0.04<=ate.post1a[,3]),mean(ate.post1a[,3]-ate.post1a[,2])),
                 c(mean((ate.post1b[,1]-0.04)^2),mean(ate.post1b[,2]<=0.04 & 0.04<=ate.post1b[,3]),mean(ate.post1b[,3]-ate.post1b[,2])),
                 c(mean((ate.post2[,1]-1)^2),mean(ate.post2[,2]<=1 & 1<=ate.post2[,3]),mean(ate.post2[,3]-ate.post2[,2])),
                 c(mean((ate.post3[,1]-ate3)^2),mean(ate.post3[,2]<=ate3 & ate3<=ate.post3[,3]),mean(ate.post3[,3]-ate.post3[,2])))
colnames(res.mat) <- c("DGP1a","DGP1b","DGP2","DGP3")
rownames(res.mat) <- c("MSE","Coverage","Int. Length")
res.mat
saveRDS(res.mat,"Results/xbcf_simulations_1.rds")
####
stopImplicitCluster()
####
