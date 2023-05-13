## Setup
## Install latest version of HighDimRD package (KR)
### devtools::install_github("kolesarm/RDHonest")
library(RDHonest)
### devtools::install_github("akreiss/HighDimRD")
library(HighDimRD)
## DGP1a
dgp <- readRDS("Data/DGP1a.rds")
s <- length(dgp) ## number of samples
for (i in 1:s)
{
    data <- dgp[[i]]
    list2env(data,globalenv())
    w1 <- fourier_basis(matrix(x),4)
    ## Estimation
    t0 <- Sys.time()
    fit  <- HighDim_rd(y,x,w1,tpc="CV" ,rd="robust")
    t1 <- Sys.time()
    write.table(t1-t0,"Results/time_kr_dgp1a.csv",append=T)
    saveRDS(fit,paste0("Results/kr_dgp1a_",i,".rds"))
}
## DGP1b
dgp <- readRDS("Data/DGP1b.rds")
s <- length(dgp) ## number of samples
for (i in 1:s)
{
    data <- dgp[[i]]
    list2env(data,globalenv())
    w1 <- fourier_basis(w,4)
    w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
    ## Estimation
    t0 <- Sys.time()
    fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
    t1 <- Sys.time()
    write.table(t1-t0,"Results/time_kr_dgp1b.csv",append=T)
    saveRDS(fit,paste0("Results/kr_dgp1b_",i,".rds"))
}
## DGP2
dgp <- readRDS("Data/DGP2.rds")
s <- length(dgp) ## number of samples
for (i in 1:s)
{
    data <- dgp[[i]]
    list2env(data,globalenv())
    w1 <- fourier_basis(w,4)
    w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
    ## Estimation
    t0 <- Sys.time()
    fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
    t1 <- Sys.time()
    write.table(t1-t0,"Results/time_kr_dgp2.csv",append=T)
    saveRDS(fit,paste0("Results/kr_dgp2_",i,".rds"))
}
## DGP1b
dgp <- readRDS("Data/DGP3.rds")
s <- length(dgp) ## number of samples
for (i in 1:s)
{
    data <- dgp[[i]]
    list2env(data,globalenv())
    w1 <- fourier_basis(w,4)
    w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
    ## Estimation
    t0 <- Sys.time()
    fit  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
    t1 <- Sys.time()
    write.table(t1-t0,"Results/time_kr_dgp3.csv",append=T)
    saveRDS(fit,paste0("Results/kr_dgp3_",i,".rds"))
}
