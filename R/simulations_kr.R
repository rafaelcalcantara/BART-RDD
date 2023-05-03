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
    ## Get interactions and expansions for KR
    w1 <- fourier_basis(w,4)
    w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
    ## Estimation
    t0 <- Sys.time()
    ate.kr  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
    t1 <- Sys.time()
    time[i] <- t1-t0
    ate.kr <- c(ate.kr$rd$Estimate[,"tau.bc"],
                ate.kr$rd$Estimate[,"tau.bc"]-1.96*ate.kr$rd$Estimate[,"se.rb"],
                ate.kr$rd$Estimate[,"tau.bc"]+1.96*ate.kr$rd$Estimate[,"se.rb"])
    results.cgs[i,] <- ate.kr
}
## KR data
results.kr <- matrix(0,s,3)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for KR data"))
    dgp(n,p,"kr")
    ## Don't get interactions and expansions for KR: already p=200 and makes it slower
    ## Estimation
    t0 <- Sys.time()
    ate.kr  <- HighDim_rd(y,x,w,tpc="CV" ,rd="robust")
    t1 <- Sys.time()
    time[s+i] <- t1-t0
    ate.kr <- c(ate.kr$rd$Estimate[,"tau.bc"],
                ate.kr$rd$Estimate[,"tau.bc"]-1.96*ate.kr$rd$Estimate[,"se.rb"],
                ate.kr$rd$Estimate[,"tau.bc"]+1.96*ate.kr$rd$Estimate[,"se.rb"])
    results.kr[i,] <- ate.kr
}
## FH data
results.fh <- matrix(0,s,3)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for FH data"))
    dgp(n,p,"fh")
    ## Get interactions and expansions for KR
    w1 <- fourier_basis(w,4)
    w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
    ## Estimation
    t0 <- Sys.time()
    ate.kr  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
    t1 <- Sys.time()
    time[2*s+i] <- t1-t0
    ate.kr <- c(ate.kr$rd$Estimate[,"tau.bc"],
                ate.kr$rd$Estimate[,"tau.bc"]-1.96*ate.kr$rd$Estimate[,"se.rb"],
                ate.kr$rd$Estimate[,"tau.bc"]+1.96*ate.kr$rd$Estimate[,"se.rb"])
    results.fh[i,] <- ate.kr
}
## Save results
results <- list(CGS=results.cgs,
                KR=results.kr,
                FH=results.fh)
saveRDS(results,"R/results_kr.rds")
saveRDS(time,"R/time_kr.rds")
