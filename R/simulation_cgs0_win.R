## Setup
set.seed(0)
library(bayesrdd)
library(doParallel)
setwd("~/../Git/BART-RDD")
if (!dir.exists("Results")) dir.create("Results") ## Create results folder
if (length(list.files("Results")[grep("cgs0_",list.files("Results"))])!=0) ## Clean up folder
{
  files <- paste0("Results/",list.files("Results")[grep("cgs0_",list.files("Results"))])
  for (i in files) file.remove(i)
}
### Parallelization
no_cores <- detectCores()-1
## Parameters
p <- c(.9,.1)
mz <- c(4, 3)
mztau <- c(3, 2)
lamstmean0_ <- rep(1,2)
lamstsd0_ <- 5*rep(3,2)
s2mean0_ <- .3
s2sd0_ <- 1
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  library(bayesrdd)
  ys <- data$y[,i]
  xs <- data$x[,i]
  bayesrdd::bayesrddest(y = ys,
                        z = xs,
                        tau = tau,
                        p = p,
                        mz = mz,
                        mztau = mztau,
                        beta0_ = rep(0,4),
                        lamstmean0_ = lamstmean0_,
                        lamstsd0_ = lamstsd0_,
                        s2mean0_ = s2mean0_,
                        s2sd0_ = s2sd0_,
                        distribution = "gaussian",
                        hetero = FALSE);
}
##
### BEGIN LOOP
files <- length(list.files("Data"))
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  n <- nrow(data$y)
  s <- ncol(data$y)
  tau <- data$c
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit)
  })
  stopCluster(cl)
  print(time)
  saveRDS(list(results=out,time=time),paste0("Results/cgs0_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")