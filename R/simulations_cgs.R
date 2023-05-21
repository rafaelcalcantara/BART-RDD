## Setup
library(parallel)
library(foreach)
library(doParallel)
## install.packages("http://apps.olin.wustl.edu/faculty/chib/rpackages/rdd/bayesrdd_1.0.zip",repo=NULL,source=T)
library(bayesrdd)
c <- 0
burn <- 100
nsamples <- 500
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## DGP1a
P = c(.70,.30)
mz = c(3, 3)
mztau = c(3, 2)
lamstmean0_ = c(1,1)
lamstsd0_ = 5*c(1,1)
s2mean0_ = .3
s2sd0_ = 1
nudgp = 3
nuval = 3
dgp <- readRDS("Data/DGP1.rds")
fit.dgp <- function(s,dgp)
{
  foreach(i=1:s,.multicombine=T,.packages = "bayesrdd",.export=c("burn","nsamples","P","mz","mztau","lamstmean0_","lamstsd0_","s2mean0_","s2sd0_","nudgp","nuval")) %do%
    {
      print(paste0("Simulation ",i," for DGP1a"))
      data <- dgp[[i]]
      list2env(data,globalenv())
      t0 <- Sys.time()
      fit  <- bayesrddest(y = y,
                          z = x,
                          tau = c,
                          p = P,
                          mz = mz,
                          mztau = mztau,
                          beta0_ = rep(0,4),
                          lamstmean0_ = lamstmean0_,
                          lamstsd0_ = lamstsd0_,
                          d = c(1,1),
                          s2mean0_ = s2mean0_,
                          s2sd0_ = s2sd0_,
                          distribution = "t",
                          nuval = nuval,
                          hetero = FALSE,
                          n0=burn,
                          m=nsamples)
      t1 <- Sys.time()
      saveRDS(list(fit=fit,time=t1-t0),paste0("Results/cgs_dgp1a_",i,".rds"))
    }
}
fit.dgp(s,dgp)
## DGP1b
P = c(.70,.30)
mz = c(3, 3)
mztau = c(3, 2)
lamstmean0_ = rep(1,12)
lamstsd0_ = 5*rep(1,12)
s2mean0_ = .3
s2sd0_ = 1
nudgp = 3
nuval = 3
fit.dgp <- function(s,dgp)
{
  foreach(i=1:s,.multicombine=T,.packages = "bayesrdd",.export=c("burn","nsamples","P","mz","mztau","lamstmean0_","lamstsd0_","s2mean0_","s2sd0_","nudgp","nuval")) %do%
    {
      print(paste0("Simulation ",i," for DGP1b"))
      data <- dgp[[i]]
      list2env(data,globalenv())
      t0 <- Sys.time()
      fit  <- bayesrddest(y = y,
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
      t1 <- Sys.time()
      saveRDS(list(fit=fit,time=t1-t0),paste0("Results/cgs_dgp1b_",i,".rds"))
    }
}
fit.dgp(s,dgp)
## DGP2
dgp <- readRDS("Data/DGP2.rds")
fit.dgp <- function(s,dgp)
{
  foreach(i=1:s,.multicombine=T,.packages = "bayesrdd",.export=c("burn","nsamples","P","mz","mztau","lamstmean0_","lamstsd0_","s2mean0_","s2sd0_","nudgp","nuval")) %do%
    {
      print(paste0("Simulation ",i," for DGP2"))
      data <- dgp[[i]]
      list2env(data,globalenv())
      t0 <- Sys.time()
      fit  <- bayesrddest(y = y,
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
      t1 <- Sys.time()
      saveRDS(list(fit=fit,time=t1-t0),paste0("Results/cgs_dgp2_",i,".rds"))
    }
}
fit.dgp(s,dgp)
## DGP3
dgp <- readRDS("Data/DGP3.rds")
fit.dgp <- function(s,dgp)
{
  foreach(i=1:s,.multicombine=T,.packages = "bayesrdd",.export=c("burn","nsamples","P","mz","mztau","lamstmean0_","lamstsd0_","s2mean0_","s2sd0_","nudgp","nuval")) %do%
    {
      print(paste0("Simulation ",i," for DGP3"))
      data <- dgp[[i]]
      list2env(data,globalenv())
      t0 <- Sys.time()
      fit  <- bayesrddest(y = y,
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
      t1 <- Sys.time()
      saveRDS(list(fit=fit,time=t1-t0),paste0("Results/cgs_dgp3_",i,".rds"))
    }
}
fit.dgp(s,dgp)
####
stopImplicitCluster()
