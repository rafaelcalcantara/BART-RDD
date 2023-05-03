set.seed(0)
## Simulations for XBCF-RDD
source("R/simulations_setup.R")
time <- rep(0,3*s)
## CGS data
results.cgs <- matrix(0,s,3)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for CGS data"))
    dgp(n,p,"cgs")
    t0 <- Sys.time()
    results.cgs[i,] <- rdd.x(y,x,w)
    t1 <- Sys.time()
    time[i] <- t1-t0
}
## KR data
results.kr <- matrix(0,s,3)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for KR data"))
    dgp(n,p,"kr")
    t0 <- Sys.time()
    results.kr[i,] <- rdd.x(y,x,w)
    t1 <- Sys.time()
    time[s+i] <- t1-t0
}
## FH data
results.fh <- matrix(0,s,3)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for FH data"))
    dgp(n,p,"fh")
    t0 <- Sys.time()
    results.fh[i,] <- rdd.x(y,x,w)
    t1 <- Sys.time()
    time[2*s+i] <- t1-t0
}
## Save results
results <- list(CGS=results.cgs,
                KR=results.kr,
                FH=results.fh)
saveRDS(results,"R/results_fh.rds")
saveRDS(time,"R/time_fh.rds")
