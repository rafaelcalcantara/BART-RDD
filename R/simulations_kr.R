## Setup
library(parallel)
library(foreach)
library(doParallel)
## Install latest version of HighDimRD package (KR)
### devtools::install_github("kolesarm/RDHonest")
library(RDHonest)
### devtools::install_github("akreiss/HighDimRD")
library(HighDimRD)
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## DGP1a
dgp <- readRDS("Data/DGP1.rds")
fit.dgp <- function(s,dgp)
{
    foreach(i=1:s,.multicombine=T) %dopar%
    {
        data <- dgp[[i]]
        list2env(data,globalenv())
        w1 <- fourier_basis(matrix(x),4)
        ## Estimation
        t0 <- Sys.time()
        fit  <- HighDim_rd(y,x,w1,tpc="CV" ,rd="robust")
        t1 <- Sys.time()
        saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_dgp1a_",i,".rds"))
    }
}
fit.dgp(s,dgp)
## DGP1b
fit.dgp <- function(s,dgp)
{
    foreach(i=1:s,.multicombine=T) %dopar%
    {
        data <- dgp[[i]]
        list2env(data,globalenv())
        w1 <- fourier_basis(matrix(x),4)
        w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
        ## Estimation
        t0 <- Sys.time()
        fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
        t1 <- Sys.time()
        saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_dgp1b_",i,".rds"))
    }
}
fit.dgp(s,dgp)
## DGP2
dgp <- readRDS("Data/DGP2.rds")
fit.dgp <- function(s,dgp)
{
    foreach(i=1:s,.multicombine=T) %dopar%
    {
        data <- dgp[[i]]
        list2env(data,globalenv())
        w1 <- fourier_basis(matrix(x),4)
        w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
        ## Estimation
        t0 <- Sys.time()
        fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
        t1 <- Sys.time()
        saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_dgp2_",i,".rds"))
    }
}
fit.dgp(s,dgp)
## DGP3
dgp <- readRDS("Data/DGP3.rds")
fit.dgp <- function(s,dgp)
{
    foreach(i=1:s,.multicombine=T) %dopar%
    {
        data <- dgp[[i]]
        list2env(data,globalenv())
        w1 <- fourier_basis(matrix(x),4)
        w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
        ## Estimation
        t0 <- Sys.time()
        fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
        t1 <- Sys.time()
        saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_dgp3_",i,".rds"))
    }
}
fit.dgp(s,dgp)
####
stopImplicitCluster()
