## ADD COVARIATES TO CGS ESTIMATION
set.seed(0)
## Simulations for XBCF-RDD
source("R/simulations_setup.R")
P = c(.70,.30);
mz = c(3, 3);
mztau = c(3, 2);
lamstmean0_ = rep(1,12);
lamstsd0_ = 5*rep(1,12);
s2mean0_ = .3;
s2sd0_ = 1;
nudgp = 3;
nuval = 3;
burn <- 100
nsamples <- 100
time <- rep(0,3*s)
## CGS data
results.cgs <- matrix(0,s,4)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for CGS data"))
    dgp(n,p,"cgs")
    t0 <- Sys.time()
    ate.cgs <- bayesrddest(y = y,
                           z = x,
                           W = w,
                           mw = rep(5,10), ## as per help file
                           tau = c,
                           p = P,
                           mz = mz,
                           mztau = mztau,
                           beta0_ = matrix(0,nrow=24),
                           lamstmean0_ = lamstmean0_,
                           lamstsd0_ = lamstsd0_,
                           d = rep(1,length(lamstsd0_)),
                           s2mean0_ = s2mean0_,
                           s2sd0_ = s2sd0_,
                           distribution = "t",
                           nuval = nuval,
                           hetero = FALSE,
                           n0=burn,
                           m=nsamples)
    ate.cgs <- c(mean(ate.cgs$atem),quantile(ate.cgs$atem,c(.025,.975)))
    results.cgs[i,1:3] <- ate.cgs
    results.cgs[i,4] <- ate
    t1 <- Sys.time()
    time[i] <- t1-t0
}
## KR data
results.kr <- matrix(0,s,4)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for KR data"))
    dgp(n,p,"kr")
    t0 <- Sys.time()
    ate.cgs <- bayesrddest(y = y,
                           z = x,
                           W = w,
                           mw = rep(5,10), ## as per help file
                           tau = c,
                           p = P,
                           mz = mz,
                           mztau = mztau,
                           beta0_ = matrix(0,nrow=24),
                           lamstmean0_ = lamstmean0_,
                           lamstsd0_ = lamstsd0_,
                           d = rep(1,length(lamstsd0_)),
                           s2mean0_ = s2mean0_,
                           s2sd0_ = s2sd0_,
                           distribution = "gaussian",
                           nuval = nuval,
                           hetero = FALSE,
                           n0=burn,
                           m=nsamples)
    ate.cgs <- c(mean(ate.cgs$atem),quantile(ate.cgs$atem,c(.025,.975)))
    results.kr[i,1:3] <- ate.cgs
    results.kr[i,4] <- ate
    t1 <- Sys.time()
    time[s+i] <- t1-t0
}
## FH data
results.fh <- matrix(0,s,4)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for FH data"))
    dgp(n,p,"fh")
    t0 <- Sys.time()
    ate.cgs <- bayesrddest(y = y,
                           z = x,
                           W = w,
                           mw = rep(5,10), ## as per help file
                           tau = c,
                           p = P,
                           mz = mz,
                           mztau = mztau,
                           beta0_ = matrix(0,nrow=24),
                           lamstmean0_ = lamstmean0_,
                           lamstsd0_ = lamstsd0_,
                           d = rep(1,length(lamstsd0_)),
                           s2mean0_ = s2mean0_,
                           s2sd0_ = s2sd0_,
                           distribution = "gaussian",
                           nuval = nuval,
                           hetero = FALSE,
                           n0=burn,
                           m=nsamples)
    ate.cgs <- c(mean(ate.cgs$atem),quantile(ate.cgs$atem,c(.025,.975)))
    results.fh[i,1:3] <- ate.cgs
    results.fh[i,4] <- ate
    t1 <- Sys.time()
    time[2*s+i] <- t1-t0
}
## Save results
results <- list(CGS=results.cgs,
                FH=results.fh,
                KR=results.kr)
saveRDS(results,"Results/results_cgs.rds")
saveRDS(time,"Results/time_cgs.rds")
