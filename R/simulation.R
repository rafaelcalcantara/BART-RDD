## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
library(rdrobust)
devtools::install_github("kolesarm/RDHonest")
library(RDHonest)
devtools::install_github("akreiss/HighDimRD")
library(HighDimRD)
s      <- 1000
sample <- c(500,1000)
model  <- 1:6
xi <- nu <- kappa <- c(0.25,2)
## Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## BART Prior settings
c             <- 0
Omin          <- 10
h             <- 0.1
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 150
burnin        <- 50
p_categorical <- 1
## BART-RDD
fit <- function(s,n,m,xi,nu,kappa)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print("BART-RDD simulations")
            print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa))
            print(paste0("Simulation ",i))
            data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",i))
            fit <- XBCF.rd(data$Y, data$W, data$X, c,
                           Owidth = h, Omin = Omin, Opct = Opct,
                           num_trees_mod = m, num_trees_con = m,
                           num_cutpoints = n,
                           num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical,
                           p_categorical_mod = p_categorical,
                           tau_con = 2*var(data$Y)/m,
                           tau_mod = 0.5*var(data$Y)/m, parallel=F)
            test <- -h<=data$X & data$X<=h
            pred <- predict.XBCFrd(fit,data$W[test,],rep(0,sum(test)))
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}
###
for (i in sample)
{
    for (j in model)
    {
        for (k in xi)
        {
            for (l in nu)
            {
                for (m in kappa)
                {
                    bart.rdd <- fit(s,i,j,k,l,m)
                    saveRDS(paste0("Results/bart_rdd_",i,"_",j,"_",k,"_",l,"_",m))
                }
            }
        }
    }
}
## CCT
fit <- function(s,n,m,xi,nu,kappa)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print("CCT simulations")
            print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa))
            print(paste0("Simulation ",i))
            data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",i))
            fit <- rdrobust(data$Y,data$X,0,covs=data$W)
            c(fit$coef[1],fit$ci[3,])
        }
}
###
for (i in sample)
{
    for (j in model)
    {
        for (k in xi)
        {
            for (l in nu)
            {
                for (m in kappa)
                {
                    cct <- fit(s,i,j,k,l,m)
                    saveRDS(paste0("Results/cct_",i,"_",j,"_",k,"_",l,"_",m))
                }
            }
        }
    }
}
## KR
fit <- function(s,n,m,xi,nu,kappa)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print("KR simulations")
            print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa))
            print(paste0("Simulation ",i))
            data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",i))
            w1 <- fourier_basis(matrix(data$X),4)
            w_HighDim <- cbind(data$W,interaction_terms(data$W),w1,interaction_terms(w1))
            ## Estimation
            fit  <- HighDim_rd(data$Y,data$X,w_HighDim,
                               tpc="CV",rd="robust")
            c(fit$rd[["Estimate"]][,"tau.bc"],fit$rd[["ci"]]["Robust",])
        }
}
###
for (i in sample)
{
    for (j in model)
    {
        for (k in xi)
        {
            for (l in nu)
            {
                for (m in kappa)
                {
                    kr <- fit(s,i,j,k,l,m)
                    saveRDS(paste0("Results/kr_",i,"_",j,"_",k,"_",l,"_",m))
                }
            }
        }
    }
}
## BART1
fit <- function(s,n,m,xi,nu,kappa)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print("BART1 simulations")
            print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa))
            print(paste0("Simulation ",i))
            data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",i))
            X <- cbind(data$X,data$W,data$Z)
            fit <- XBART(data$Y, X, num_trees = m,
                         num_cutpoints = n, num_sweeps = num_sweeps,
                         burnin = burnin, Nmin = Nmin,
                         p_categorical_con = p_categorical+1,
                         ## add 1 for Z
                         p_categorical_mod = p_categorical+1,
                         tau_con = 2*var(y)/m,
                         tau_mod = 0.5*var(y)/m, parallel=F)
            pred1 <- predict.XBART(fit,cbind(rep(0,n),data$W,rep(1,n)))[,(burnin+1):num_sweeps]
            pred0 <- predict.XBART(fit,cbind(rep(0,n),data$W,rep(0,n)))[,(burnin+1):num_sweeps]
            pred1-pred0
        }
}
###
for (i in sample)
{
    for (j in model)
    {
        for (k in xi)
        {
            for (l in nu)
            {
                for (m in kappa)
                {
                    bart1 <- fit(s,i,j,k,l,m)
                    saveRDS(paste0("Results/bart1_",i,"_",j,"_",k,"_",l,"_",m))
                }
            }
        }
    }
}
## BART2
fit <- function(s,n,m,xi,nu,kappa)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print("BART2 simulations")
            print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa))
            print(paste0("Simulation ",i))
            data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",i))
            X <- cbind(data$X,data$W)
            fit0 <- XBART(data$Y[!data$Z], X[!data$Z,],
                          num_trees = m,
                          num_cutpoints = n,
                          num_sweeps = num_sweeps,
                          burnin = burnin, Nmin = Nmin,
                          p_categorical_con = p_categorical,
                          p_categorical_mod = p_categorical,
                          tau_con = 2*var(data$Y[!data$Z])/m,
                          tau_mod = 0.5*var(data$Y[!data$Z])/m,
                          parallel=F)
            fit1 <- XBART(data$Y[data$Z], X[data$Z,],
                          num_trees = m,
                          num_cutpoints = n,
                          num_sweeps = num_sweeps,
                          burnin = burnin, Nmin = Nmin,
                          p_categorical_con = p_categorical,
                          p_categorical_mod = p_categorical,
                          tau_con = 2*var(data$Y[data$Z])/m,
                          tau_mod = 0.5*var(data$Y[data$Z])/m,
                          parallel=F)
            pred1 <- predict.XBART(fit1,cbind(rep(0,sum(data$Z)),data$W))
            pred0 <- predict.XBART(fit0,cbind(rep(0,sum(!data$Z)),data$W))
            colMeans(pred1[,(burnin+1):num_sweeps])-colMeans(pred0[,(burnin+1):num_sweeps])
        }
}
###
for (i in sample)
{
    for (j in model)
    {
        for (k in xi)
        {
            for (l in nu)
            {
                for (m in kappa)
                {
                    bart1 <- fit(s,i,j,k,l,m)
                    saveRDS(paste0("Results/bart2_",i,"_",j,"_",k,"_",l,"_",m))
                }
            }
        }
    }
}
BCF
fit <- function(s,n,m,xi,nu,kappa)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print("BCF simulations")
            print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa))
            print(paste0("Simulation ",i))
            data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",i))
            X <- cbind(data$X,data$W)
            fit <- XBCF.discrete(data$Y, data$Z, X_con=X, X_mod=X,
                                 num_trees_con = m,
                                 num_trees_mod = m,
                                 num_cutpoints = n,
                                 num_sweeps = num_sweeps,
                                 burnin = burnin, Nmin = Nmin,
                                 p_categorical_con = p_categorical,
                                 p_categorical_mod = p_categorical,
                                 tau_con = 2*var(data$Y)/m,
                                 tau_mod = 0.5*var(data$Y)/m,
                                 parallel=F)
            pred <- predict.XBCFdiscrete(fit,cbind(rep(0,n),data$W),
                                         cbind(rep(0,n),data$W),
                                         data$Z,pihat=data$Z,
                                         burnin=burnin)
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}

for (i in sample)
{
    for (j in model)
    {
        for (k in xi)
        {
            for (l in nu)
            {
                for (m in kappa)
                {
                    bart1 <- fit(s,i,j,k,l,m)
                    saveRDS(paste0("Results/bcf_",i,"_",j,"_",k,"_",l,"_",m))
                }
            }
        }
    }
}
