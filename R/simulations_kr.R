## Setup
library(parallel)
library(foreach)
library(doParallel)
## Install latest version of HighDimRD package (KR)
devtools::install_github("kolesarm/RDHonest")
library(RDHonest)
devtools::install_github("akreiss/HighDimRD")
library(HighDimRD)
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## DGP1
for (p in c(4,6,10))
{
    ## DGP1a
    dgp <- readRDS(paste0("Data/DGP1_",p,".rds"))
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 1a, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w1,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                saveRDS(list(fit=fit,time=dt),paste0("Results/kr_1a_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,p,dgp)
    ## DGP1b
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 1b, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- d
                saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_1b_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,dgp)
}
## DGP2
for (p in c(4,6,10))
{
    ## DGP2a
    dgp <- readRDS(paste0("Data/DGP2_",p,".rds"))
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 2a, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w1,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                saveRDS(list(fit=fit,time=dt),paste0("Results/kr_2a_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,p,dgp)
    ## DGP2b
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 2b, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- d
                saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_2b_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,dgp)
}
## DGP3
for (p in c(4,6,10))
{
    ## DGP3a
    dgp <- readRDS(paste0("Data/DGP3_",p,".rds"))
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 3a, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w1,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                saveRDS(list(fit=fit,time=dt),paste0("Results/kr_3a_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,p,dgp)
    ## DGP3b
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 3b, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- d
                saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_3b_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,dgp)
}
## DGP4
for (p in c(4,6,10))
{
    ## DGP4a
    dgp <- readRDS(paste0("Data/DGP4_",p,".rds"))
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 4a, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w1,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                saveRDS(list(fit=fit,time=dt),paste0("Results/kr_4a_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,p,dgp)
    ## DGP4b
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 4b, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- d
                saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_4b_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,dgp)
}
## DGP5
for (p in c(4,6,10))
{
    ## DGP5a
    dgp <- readRDS(paste0("Data/DGP5_",p,".rds"))
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 5a, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w1,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                saveRDS(list(fit=fit,time=dt),paste0("Results/kr_5a_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,p,dgp)
    ## DGP5b
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 5b, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- d
                saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_5b_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,dgp)
}
## DGP6
for (p in c(4,6,10))
{
    ## DGP6a
    dgp <- readRDS(paste0("Data/DGP6_",p,".rds"))
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 6a, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w1,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                saveRDS(list(fit=fit,time=dt),paste0("Results/kr_6a_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,p,dgp)
    ## DGP6b
    fit.dgp <- function(s,p,dgp)
    {
        foreach(i=1:s,.multicombine=T) %dopar%
            {
                print(paste0("Simulation ",i," for DGP 6b, ",p," covariates"))
                data <- dgp[[i]]
                list2env(data,globalenv())
                w1 <- fourier_basis(matrix(x),4)
                w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
                ## Estimation
                t0 <- Sys.time()
                fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
                t1 <- Sys.time()
                dt <- d
                saveRDS(list(fit=fit,time=t1-t0),paste0("Results/kr_6b_",p,"_",i,".rds"))
            }
    }
    fit.dgp(s,dgp)
}
####
stopImplicitCluster()
