## Setup
library(parallel)
library(foreach)
library(doParallel)
library(rdrobust)
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## Functions
ckt.fit <- function(s,y,x,w)
{
    foreach(i=1:s,.multicombine=T,.export="c") %dopar%
        {
            rdrobust(y[,i],x[,i],c=c,covs=w[[i]])
        }
}
## DGP1
list2env(readRDS("Data/dgp1.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
ckt.time <- system.time({
    ckt <- ckt.fit(s,y,x,w)
})
saveRDS(list(fit=ckt,time=ckt.time),"Results/ckt_1.rds")
rm(y,x,w,z,s,n,ckt.time,ckt)
## DGP2
list2env(readRDS("Data/dgp2.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
ckt.time <- system.time({
    ckt <- ckt.fit(s,y,x,w)
})
saveRDS(list(fit=ckt,time=ckt.time),"Results/ckt_2.rds")
rm(y,x,w,z,s,n,ckt.time,ckt)
## DGP3
list2env(readRDS("Data/dgp3.rds"),globalenv())
s <- ncol(x)
n <- nrow(x)
ckt.time <- system.time({
    ckt <- ckt.fit(s,y,x,w)
})
saveRDS(list(fit=ckt,time=ckt.time),"Results/ckt_3.rds")
rm(y,x,w,z,s,n,ckt.time,ckt)
##
print("Done!")
