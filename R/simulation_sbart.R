## Setup
set.seed(0)
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
### Parameters
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 3
### Functions
fit <- function(s,y,x,w,z)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("Sample: ",i))
            ys <- y[,i]
            ws <- w[[i]]
            xs <- x[,i]
            zs <- z[,i]
            fit <- XBART(ys, cbind(xs,ws,zs), num_trees = ntrees,
                         num_cutpoints = n, num_sweeps = num_sweeps,
                         burnin = burnin, Nmin = Nmin,
                         p_categorical = p_categorical,
                         tau = var(ys)/ntrees, parallel=F)
            pred1 <- predict.XBART(fit,cbind(rep(0,n),ws,rep(1,n)))[,(burnin+1):num_sweeps]
            pred0 <- predict.XBART(fit,cbind(rep(0,n),ws,rep(0,n)))[,(burnin+1):num_sweeps]
            pred1-pred0
        }
}
## DGP1
list2env(readRDS("Data/dgp1.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
bart.time1 <- system.time({
    bart.fit1 <- fit(s,y,x,w,z)
})
saveRDS(list(fit=bart.fit1,time=bart.time1),"Results/sbart_1.rds")
rm(s,n,x,y,w,z,bart.time1,bart.fit1)
## DGP2
list2env(readRDS("Data/dgp2.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
bart.time2 <- system.time({
    bart.fit2 <- fit(s,y,x,w,z)
})
saveRDS(list(fit=bart.fit2,time=bart.time2),"Results/sbart_2.rds")
rm(s,n,x,y,w,z,bart.time2,bart.fit2)
## DGP3
list2env(readRDS("Data/dgp3.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
bart.time3 <- system.time({
    bart.fit3 <- fit(s,y,x,w,z)
})
saveRDS(list(fit=bart.fit3,time=bart.time3),"Results/sbart_3.rds")
rm(s,n,x,y,w,z,bart.time3,bart.fit3)
##
print("Done!")
