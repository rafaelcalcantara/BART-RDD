## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
###
fit.xbcf <- function(y,w,x)
{
    t0 <- Sys.time()
    h <- 0.15
    fit <- XBCF.rd(y, W=w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
                   num_trees_mod = m, num_trees_con = m,
                   num_cutpoints = n, num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   tau_con = 2*var(y)/m,
                   tau_mod = 0.5*var(y)/m, parallel=F)
    test <- -h<=x & x<=h
    pred <- predict.XBCFrd(fit,w[test,],rep(0,sum(test)))
    post <- pred$tau.adj[,(burnin+1):num_sweeps]
    ## ind.post <- predict.XBCFrd(fit,w,x)
    ## ind.post <- ind.post$tau.adj[,(burnin+1):num_sweeps]
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    return(list(ate.post=post,Owidth=h,time=dt))
}
fit.1a <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 1a, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,NULL,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_1a_",p,"_",i,".rds"))
        }
}
fit.1b <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 1b, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_1b_",p,"_",i,".rds"))
        }
}
fit.2a <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 2a, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,NULL,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_2a_",p,"_",i,".rds"))
        }
}
fit.2b <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 2b, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_2b_",p,"_",i,".rds"))
        }
}
fit.3a <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 3a, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,NULL,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_3a_",p,"_",i,".rds"))
        }
}
fit.3b <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 3b, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_3b_",p,"_",i,".rds"))
        }
}
fit.4a <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 4a, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,NULL,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_4a_",p,"_",i,".rds"))
        }
}
fit.4b <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 4b, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_4b_",p,"_",i,".rds"))
        }
}
fit.5a <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 5a, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,NULL,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_5a_",p,"_",i,".rds"))
        }
}
fit.5b <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 5b, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_5b_",p,"_",i,".rds"))
        }
}
fit.6a <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 6a, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,NULL,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_6a_",p,"_",i,".rds"))
        }
}
fit.6b <- function(s,p,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 6b, ",p," covariates"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_6b_",p,"_",i,".rds"))
        }
}
fit.7 <- function(s,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("XBCF: Simulation ",i," for DGP 7"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_7_",i,".rds"))
        }
}
####
c             <- 0
Omin          <- 10
## Omin          <- as.integer(0.03*n)
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 0
num_cutpoints <- n
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## DGP1
for (p in c(4,6,10))
{
    dgp <- readRDS(paste0("Data/DGP1_",p,".rds"))
    fit.1a(s,p,dgp)
    ## fit.1b(s,p,dgp)
}
## DGP2
for (p in c(4,6,10))
{
    dgp <- readRDS(paste0("Data/DGP2_",p,".rds"))
    fit.2a(s,p,dgp)
    ## fit.2b(s,p,dgp)
}
## DGP3
for (p in c(4,6,10))
{
    dgp <- readRDS(paste0("Data/DGP3_",p,".rds"))
    fit.3a(s,p,dgp)
    ## fit.3b(s,p,dgp)
}
## DGP4
for (p in c(4,6,10))
{
    dgp <- readRDS(paste0("Data/DGP4_",p,".rds"))
    fit.4a(s,p,dgp)
    ## fit.4b(s,p,dgp)
}
## DGP5
for (p in c(4,6,10))
{
    dgp <- readRDS(paste0("Data/DGP5_",p,".rds"))
    fit.5a(s,p,dgp)
    ## fit.5b(s,p,dgp)
}
## DGP6
for (p in c(4,6,10))
{
    dgp <- readRDS(paste0("Data/DGP6_",p,".rds"))
    fit.6a(s,p,dgp)
    ## fit.6b(s,p,dgp)
}
## DGP7
## dgp <- readRDS("Data/DGP7.rds")
## fit.7(s,dgp)
####
stopImplicitCluster()
