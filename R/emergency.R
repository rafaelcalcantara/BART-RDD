rm(list=ls())
n <- 500
p <- 10
c <- 0
s <- 1000 # samples
setwd("~/Documents/Git/XBCF-RDD/")
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
results.fh.xbcf <- matrix(0,s,4)
for (i in 1:s)
{
    print(paste0("Simulation ",i," for FH data"))
    dgp(n,p,"fh")
    t0 <- Sys.time()
    Owidth <- findOwidth(0.01)
    results.fh.xbcf[i,1:3] <- pred.ate.xbcf()
    results.fh.xbcf[i,4] <- ate
    t1 <- Sys.time()
    time[3*s+i] <- t1-t0
}
##
results.fh.kr <- matrix(0,s,4)
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
    time[3*s+i] <- t1-t0
    ate.kr <- c(ate.kr$rd$Estimate[,"tau.bc"],
                ate.kr$rd$Estimate[,"tau.bc"]-1.96*ate.kr$rd$Estimate[,"se.rb"],
                ate.kr$rd$Estimate[,"tau.bc"]+1.96*ate.kr$rd$Estimate[,"se.rb"])
    results.fh.kr[i,1:3] <- ate.kr
    results.fh.kr[i,4] <- ate
}
saveRDS(list(xbcf=results.fh.xbcf,kr=results.fh.kr),
        "Results/results_fh.rds")
