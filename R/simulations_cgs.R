## Setup
library(parallel)
library(foreach)
library(doParallel)
## install.packages("http://apps.olin.wustl.edu/faculty/chib/rpackages/rdd/bayesrdd_1.0.zip",repo=NULL,source=T)
library(bayesrdd)
c <- 0
burn <- 100
nsamples <- 150
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## j = 1; i = 4; k = 1; ab = "a"; p = 4
###
for (j in 1:6)
{
  for (i in c(4,6,10))
  {
    dgp <- readRDS(paste0("Data/DGP",j,"_",i,".rds"))
    fit.dgp <- function(s,p,dgp,ab)
    {
      foreach(k=1:s,.multicombine=T,.packages = "bayesrdd",.export=c("burn","nsamples")) %do%
        {
          P = c(.70,.30)
          mz = c(5,5)
          mztau = c(5,5)
          s2mean0_ = .3
          s2sd0_ = 1
          print(paste0("Simulation ",k," for DGP",j,ab,", ",p," covariates"))
          data <- dgp[[k]]
          list2env(data,globalenv())
          t0 <- Sys.time()
          if (ab=="a")
          {
            lamstmean0_ = rep(1,2)
            lamstsd0_ = 5*rep(1,2)
            fit  <- bayesrddest(y = y,
                                z = x,
                                tau = c,
                                p = P,
                                mz = mz,
                                mztau = mztau,
                                beta0_ = rep(0,4),
                                lamstmean0_ = lamstmean0_,
                                lamstsd0_ = lamstsd0_,
                                d = rep(1,2),
                                s2mean0_ = s2mean0_,
                                s2sd0_ = s2sd0_,
                                distribution = "gaussian",
                                nuval=3,
                                hetero = FALSE,
                                n0=burn,
                                m=nsamples)
          } else
          {
            lamstmean0_ = rep(1,p+2)
            lamstsd0_ = 5*rep(1,p+2)
            fit  <- bayesrddest(y = y,
                                z = x,
                                W = w,
                                mw = rep(5,p), ## as per help file
                                tau = c,
                                p = P,
                                mz = mz,
                                mztau = mztau,
                                beta0_ = rep(0,4+2*p),
                                lamstmean0_ = lamstmean0_,
                                lamstsd0_ = lamstsd0_,
                                d = rep(1,p+2),
                                s2mean0_ = s2mean0_,
                                s2sd0_ = s2sd0_,
                                distribution = "gaussian",
                                nuval=3,
                                hetero = FALSE,
                                n0=burn,
                                m=nsamples)
          }
          t1 <- Sys.time()
          dt <- difftime(t1,t0)
          saveRDS(list(ate.post=fit$atem,time=dt),paste0("Results/cgs",j,"_",ab,"_",p,"_",k,".rds"))
        }
    }
    fit.dgp(s,i,dgp,"a")
    fit.dgp(s,i,dgp,"b")
  }
}
####
stopImplicitCluster()
