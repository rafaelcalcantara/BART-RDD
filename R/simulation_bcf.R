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
p_categorical <- 2
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
            fit <- XBCF.continuous(y=ys, Z=zs, X_con = cbind(xs,ws), X_mod = cbind(xs,ws),
                           num_trees_mod = ntrees,
                           num_trees_con = ntrees,
                           num_cutpoints = n,
                           num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical,
                           p_categorical_mod = p_categorical,
                           tau_con = 2*var(ys)/ntrees,
                           tau_mod = 0.5*var(ys)/ntrees)
            pred <- predict.XBCFcontinuous(fit,X_con = cbind(xs,ws), X_mod = cbind(xs,ws),Z=zs)
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}
## DGP1
list2env(readRDS("Data/dgp1.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
bart.time1 <- system.time({
    bart.fit1 <- fit(s,y,x,w,z)
})
saveRDS(list(fit=bart.fit1,time=bart.time1),"Results/bcf_1.rds")
rm(s,n,x,y,w,z,bart.time1,bart.fit1)
## DGP2
list2env(readRDS("Data/dgp2.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
bart.time2 <- system.time({
    bart.fit2 <- fit(s,y,x,w,z)
})
saveRDS(list(fit=bart.fit2,time=bart.time2),"Results/bcf_2.rds")
rm(s,n,x,y,w,z,bart.time2,bart.fit2)
## DGP3
list2env(readRDS("Data/dgp3.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
bart.time3 <- system.time({
    bart.fit3 <- fit(s,y,x,w,z)
})
saveRDS(list(fit=bart.fit3,time=bart.time3),"Results/bcf_3.rds")
rm(s,n,x,y,w,z,bart.time3,bart.fit3)
##
print("Done!")
