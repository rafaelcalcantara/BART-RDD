set.seed(0)
## Simulations for XBCF-RDD
source("R/simulations_setup.R")
time <- rep(0,4*s)
## Settings for XBCF-RDD
### Owidth: overlap bandwidth (defined per data)
### Omin: minimum number of observations inside overlap region for each leaf node
### Opct: If nb is the number of obs in node b, Opct*nb of them have to be inside overlap region
Omin          <- 2
Opct          <- 0.9
num_trees_mod <- 10
num_trees_con <- 10
num_cutpoints <- n
Nmin          <- 10
num_sweeps    <- 50
burnin        <- 20
p_categorical <- 0
## CGS data
### No covariates
results.cgs0 <- matrix(0,s,4)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for CGS data - No covariates"))
    dgp(n,p,"cgs")
    t0 <- Sys.time()
    Owidth <- findOwidth(0.01)
    results.cgs0[i,1:3] <- pred.ate.xbcf0()
    results.cgs0[i,4] <- ate
    t1 <- Sys.time()
    time[i] <- t1-t0
}
### With covariates
results.cgs <- matrix(0,s,4)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for CGS data - With covariates"))
    dgp(n,p,"cgs")
    t0 <- Sys.time()
    Owidth <- findOwidth(0.01)
    results.cgs[i,1:3] <- pred.ate.xbcf()
    results.cgs[i,4] <- ate
    t1 <- Sys.time()
    time[s+i] <- t1-t0
}
## KR data
results.kr <- matrix(0,s,4)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for KR data"))
    dgp(n,p,"kr")
    t0 <- Sys.time()
    Owidth <- findOwidth(0.01)
    results.kr[i,1:3] <- pred.ate.xbcf()
    results.kr[i,4] <- ate
    t1 <- Sys.time()
    time[2*s+i] <- t1-t0
}
## FH data
results.fh <- matrix(0,s,4)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for FH data"))
    dgp(n,p,"fh")
    t0 <- Sys.time()
    Owidth <- findOwidth(0.01)
    results.fh[i,1:3] <- pred.ate.xbcf()
    results.fh[i,4] <- ate
    t1 <- Sys.time()
    time[3*s+i] <- t1-t0
}
## Save results
results <- list(CGS0 = results.cgs0,
                CGS=results.cgs,
                FH=results.fh,
                KR=results.kr)
saveRDS(results,"Results/results_xbcf.rds")
saveRDS(time,"Results/time_xbcf.rds")
