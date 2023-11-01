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
Omin          <- 5
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 2
### Functions
opt.h.int <- function(s,h,ys,x,w,z)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            y <- ys[,i]
            print(paste0("Sample: ",i))
            fit <- XBCF.rd(y, w, x, c,
                           Owidth = h, Omin = Omin, Opct = Opct,
                           num_trees_mod = ntrees,
                           num_trees_con = ntrees,
                           num_cutpoints = n,
                           num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical,
                           p_categorical_mod = p_categorical,
                           tau_con = 2*var(y)/ntrees,
                           tau_mod = 0.5*var(y)/ntrees)
            test <- -h+c<=x & x<=h+c
            pred <- predict.XBCFrd(fit,w[test,],rep(c,sum(test)))
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}
opt.h <- function(s,y,x,w,z,h)
{
    out <- vector("list",length(h))
    for (i in 1:length(h))
    {
        print(paste0("h: ",i))
        temp <- opt.h.int(s,h[i],y,x,w,z)
        out[[i]] <- temp
    }
    return(out)
}
## h <- round(rbind(0.5*apply(x,2,sd),apply(x,2,sd),2*apply(x,2,sd)),2)
## h1 <- vector("list",s2)
## bart.prior.time <- system.time({
##     for (j in 1:s2) h1[[j]] <- opt.h(s1,ys[[j]],x[,j],w[[j]],z[,j],h[,j])
## })
## saveRDS(h1,"Results/prior_h.rds")
## h1 <- readRDS("Results/prior_h.rds")
## out <- matrix(0,s2,nrow(h))
## for (i in 1:s2)
## {
##     for (j in 1:nrow(h))
##     {
##         out[i,j] <- sqrt(mean(sapply(h1[[i]][[j]],function(x) (mean(colMeans(x))-mean(tau.prior(0,w[[i]])))^2)))
##     }
## }
## h.opt <- mapply(function(i,j) h[j,i],i=1:s2,j=apply(out,1,function(x) which(x==min(x))))
## ##
fit <- function(s,h,y,x,w,z)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("Sample: ",i))
            ys <- y[,i]
            hs <- h[i]
            ws <- w[[i]]
            xs <- x[,i]
            fit <- XBCF.rd(ys, ws, xs, c,
                           Owidth = hs, Omin = Omin, Opct = Opct,
                           num_trees_mod = ntrees,
                           num_trees_con = ntrees,
                           num_cutpoints = n,
                           num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical,
                           p_categorical_mod = p_categorical,
                           tau_con = 2*var(ys)/ntrees,
                           tau_mod = 0.5*var(ys)/ntrees)
            test <- -hs+c<=xs & xs<=hs+c
            pred <- predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}
## DGP1
list2env(readRDS("Data/dgp1.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
h <- apply(x,2,sd)*0.5
bart.time1 <- system.time({
    bart.fit1 <- fit(s,h,y,x,w,z)
})
saveRDS(list(fit=bart.fit1,time=bart.time1),"Results/bart_rdd_1.rds")
rm(s,n,x,y,w,z,h,bart.time1,bart.fit1)
## DGP2
list2env(readRDS("Data/dgp2.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
h <- apply(x,2,sd)*0.5
bart.time2 <- system.time({
    bart.fit2 <- fit(s,h,y,x,w,z)
})
saveRDS(list(fit=bart.fit2,time=bart.time2),"Results/bart_rdd_2.rds")
rm(s,n,x,y,w,z,h,bart.time2,bart.fit2)
## DGP3
list2env(readRDS("Data/dgp3.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
h <- apply(x,2,sd)*0.5
bart.time3 <- system.time({
    bart.fit3 <- fit(s,h,y,x,w,z)
})
saveRDS(list(fit=bart.fit3,time=bart.time3),"Results/bart_rdd_3.rds")
rm(s,n,x,y,w,z,h,bart.time3,bart.fit3)
##
print("Done!")
